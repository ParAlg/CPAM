#pragma once

#include <cpam/cpam.h>
#include <pam/pam.h>
// using namespace cpam;
#include <pam/get_time.h>
#include <pam/parse_command_line.h>
#include "../graphs/aspen/aspen.h"
#include "types.h"
#include "util/NSGDist.h"
#include "util/counter.h"
#include "util/distance.h"

bool stats = false;

template <typename T>
struct knn_index {
  // Pointers to the coordinates for each point.
  parlay::sequence<Tvec_point<T>*>& v;

  size_t maxDeg;
  size_t beamSize;
  double alpha;
  size_t d;

  atomic_sum_counter<size_t> total_visited;

  struct empty_weight {};
  using Graph = aspen::symmetric_graph<empty_weight>;
  using edge_tree = typename Graph::edge_tree;
  using edge_node = typename Graph::edge_node;
  using vertex_tree = typename Graph::vertex_tree;
  using vertex_node = typename Graph::vertex_node;

  Graph G;

  using tvec_point = Tvec_point<T>;
  using fvec_point = Tvec_point<float>;

  tvec_point* medoid;

  using pid = std::pair<node_id, float>;
  using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));

  knn_index(parlay::sequence<Tvec_point<T>*>& v, size_t maxDeg, size_t beamSize,
            double Alpha, size_t dim)
      : v(v), maxDeg(maxDeg), beamSize(beamSize), alpha(Alpha), d(dim) {
    std::cout << "Initialized knn_index with maxDeg = " << maxDeg
              << " beamSize = " << beamSize << " alpha = " << alpha
              << " dim = " << dim << std::endl;

    total_visited.reset();
  }

  void print_graph_status() {
    std::cout << "G.num_vertices = " << G.num_vertices()
              << " num_edges = " << G.num_edges() << std::endl;
  }

  void build_index(parlay::sequence<node_id> inserts) {
    // Find the medoid, which each beamSearch will begin from.
    node_id medoid_id = find_approx_medoid();
    G.insert_vertex_inplace(medoid_id, nullptr);
    batch_insert(inserts, 2, .02, false);
    std::cout << "G.num_vertices = " << G.num_vertices()
              << " num_edges = " << G.num_edges() << std::endl;
#ifdef STATS
    std::cout << "Total vertices visited: " << total_visited.get_value()
              << std::endl;
    std::cout << "Total distance calls: " << distance_calls.get_value()
              << std::endl;
#endif

    // compute "self-recall", i.e. recall from base points.
    // Self-recall and query recall being different is an indicator
    // of distribution shift.

    if (stats) {
      auto print_fn = [&](node_id u, node_id v, empty_weight wgh) -> bool {
        std::cout << "Neighbor = " << v << std::endl;
        return true;
      };
      auto m = G.get_vertex(medoid->id);
      std::cout << "Medoid's neighbors: " << std::endl;
      m.out_neighbors().foreach_cond(print_fn);

      compute_self_recall();
    }
  }

  void insert(node_id p) {
    parlay::sequence<node_id> inserts;
    inserts.push_back(p);
    batch_insert(inserts, 2, .02, false);
  }

  void insert(parlay::sequence<node_id> inserts, bool random_order = false) {
    batch_insert(inserts, 2, .02, random_order);
  }

  parlay::sequence<node_id> query(T* query_coords, int k, int beamSizeQ) {
    if ((k + 1) > beamSizeQ) {
      std::cout << "Error: beam search parameter Q = " << beamSizeQ
                << " same size or smaller than k = " << k << std::endl;
      abort();
    }
    auto pairs = beam_search(query_coords, beamSizeQ);
    auto& beamElts = pairs.first;
    parlay::sequence<node_id> neighbors(k);
    // Ignoring reporting the point itself for now.
    for (int j = 0; j < k; j++) {
      neighbors[j] = beamElts[j].first;
    }
    return neighbors;
  }

  void lazy_delete(parlay::sequence<node_id> deletes) {
    for (node_id p : deletes) {
      if (p < 0 || p > (node_id)v.size()) {
        std::cout << "ERROR: invalid point " << p << " given to lazy_delete"
                  << std::endl;
        abort();
      }
      if (p != medoid->id)
        delete_set.insert(p);
      else
        std::cout << "Deleting medoid not permitted; continuing" << std::endl;
    }
  }

  void lazy_delete(node_id p) {
    if (p < 0 || p > (node_id)v.size()) {
      std::cout << "ERROR: invalid point " << p << " given to lazy_delete"
                << std::endl;
      abort();
    }
    if (p == (node_id)medoid->id) {
      std::cout << "Deleting medoid not permitted; continuing" << std::endl;
      return;
    }
    delete_set.insert(p);
  }

  void consolidate_deletes() {
    // freeze the delete set before beginning a consolidate
    // later will add versioning/locks
    parlay::sequence<node_id> delete_vec;
    for (auto d : delete_set) delete_vec.push_back(d);
    std::set<unsigned> old_delete_set;
    delete_set.swap(old_delete_set);
    // auto to_consolidate = parlay::tabulate(v.size(), [&] (node_id i) {return
    // i;});
    // consolidate_deletes_internal(old_delete_set, to_consolidate);
    consolidate_deletes_simple(old_delete_set);
    remove_deleted_vertices(delete_vec);
    check_deletes_correct(old_delete_set);
  }

 private:
  std::set<node_id> delete_set;
  // p_coords: query vector coordinates
  // v: database of vectors
  // medoid: "root" of the proximity graph
  // beamSize: (similar to ef)
  // d: dimensionality of the indexed vectors

  void consolidate_deletes_simple(std::set<node_id> old_delete_set) {
    auto consolidated_vertices =
        parlay::sequence<std::tuple<node_id, edge_node*>>(v.size());

    parlay::parallel_for(0, v.size(), [&](size_t i) {
      if (old_delete_set.find(i) == old_delete_set.end()) {
        auto current_vtx = G.get_vertex(i);
        parlay::sequence<node_id> candidates;
        auto g = [&](node_id a) {
          return old_delete_set.find(a) == old_delete_set.end();
        };
        auto f = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v)) candidates.push_back(v);
          return true;
        };
        current_vtx.out_neighbors().foreach_cond(f);
        auto begin = (std::tuple<node_id, empty_weight>*)candidates.begin();
        auto tree = edge_tree(begin, begin + candidates.size());
        consolidated_vertices[i] = {i, tree.root};
        tree.root = nullptr;
      }
    });
    G.insert_vertices_batch(consolidated_vertices.size(),
                            consolidated_vertices.begin());
  }

  void consolidate_deletes_internal(std::set<node_id> old_delete_set,
                                    parlay::sequence<node_id>& to_consolidate) {
    auto consolidated_vertices =
        parlay::sequence<std::tuple<node_id, edge_node*>>(
            to_consolidate.size());
    std::vector<bool> needs_consolidate(to_consolidate.size(), true);

    parlay::parallel_for(0, to_consolidate.size(), [&](size_t i) {
      node_id index = to_consolidate[i];
      if (old_delete_set.find(index) == old_delete_set.end()) {
        // remove deleted vertices and add their out_nbh
        bool change = false;
        auto current_vtx = G.get_vertex(index);
        parlay::sequence<node_id> candidates;
        auto g = [&](node_id a) {
          return old_delete_set.find(a) == old_delete_set.end();
        };

        auto f = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v)) candidates.push_back(v);
          return true;
        };

        auto h = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v))
            candidates.push_back(v);
          else {
            change = true;
            auto vtx = G.get_vertex(v);
            vtx.out_neighbors().foreach_cond(f);
          }
          return true;
        };
        current_vtx.out_neighbors().foreach_cond(h);
        // for(auto cand : candidates){
        //   if(old_delete_set.find(cand) != old_delete_set.end()){
        //     std::cout << "ERROR: after assembling candidate list, " <<
        //     std::endl;
        //     std::cout << "vertex " << index << " candidate list contains
        //     deleted neighbor "
        //       << cand << std::endl;
        //   }
        // }
        if (change) {
          needs_consolidate[i] = true;
          if (candidates.size() <= maxDeg) {
            auto begin = (std::tuple<node_id, empty_weight>*)candidates.begin();
            auto tree = edge_tree(begin, begin + candidates.size());
            consolidated_vertices[i] = {index, tree.root};
            tree.root = nullptr;
          } else {
            parlay::sequence<node_id> new_out_2(
                maxDeg, std::numeric_limits<node_id>::max());
            auto output_slice = parlay::make_slice(new_out_2.begin(),
                                                   new_out_2.begin() + maxDeg);
            robustPrune(v[index], index, candidates, alpha, output_slice,
                        false);
            size_t deg = size_of(output_slice);
            // for(size_t j=0; j<deg; j++){
            //   if(old_delete_set.find(new_out_2[j]) != old_delete_set.end()){
            //     std::cout << "ERROR: after pruning candidate list, " <<
            //     std::endl;
            //     std::cout << "vertex " << index << " candidate list contains
            //     deleted neighbor "
            //       << new_out_2[j] << std::endl;
            //   }
            // }
            auto begin = (std::tuple<node_id, empty_weight>*)new_out_2.begin();
            auto tree = edge_tree(begin, begin + deg);
            consolidated_vertices[i] = {index, tree.root};
            tree.root = nullptr;
          }
        } else {
          needs_consolidate[i] = false;
        }
      }
    });
    // auto vtx = G.get_vertex(1348);
    // auto f = [&] (node_id v, node_id u, empty_weight wgh){
    //   std::cout << u << std::endl;
    //   return true;
    // };
    // std::cout << "before" << std::endl;
    // vtx.out_neighbors().foreach_cond(f);
    auto filtered_vertices =
        parlay::pack(consolidated_vertices, needs_consolidate);
    // std::cout << filtered_vertices.size() << std::endl;
    // for(auto fv : filtered_vertices){if(get<0>(fv) == 1348) std::cout <<
    // "HERE" << std::endl;}
    // std::cout << std::endl;
    // std::cout << "after" << std::endl;
    // vtx.out_neighbors().foreach_cond(f);
    G.insert_vertices_batch(filtered_vertices.size(),
                            filtered_vertices.begin());
  }

  void remove_deleted_vertices(parlay::sequence<node_id>& delete_vec) {
    G.delete_vertices_batch(delete_vec.size(), delete_vec.begin());
  }

  void check_deletes_correct(std::set<node_id>& old_delete_set) {
    parlay::parallel_for(0, v.size(), [&](size_t i) {
      if (old_delete_set.find(i) == old_delete_set.end()) {
        auto current_vtx = G.get_vertex(i);
        auto g = [&](node_id a) {
          return (old_delete_set.find(a) != old_delete_set.end());
        };
        parlay::sequence<node_id> remaining_nbh;
        auto f = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v)) {
            std::cout << "ERROR: vertex " << u << " has deleted neighbor " << v
                      << std::endl;
          }
          return true;
        };
        current_vtx.out_neighbors().foreach_cond(f);
      }
    });
  }

  std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
      T* p_coords, int beamSize) {
    // initialize data structures
    // pid = std::pair<node_id, float>
    std::vector<pid> visited;
    parlay::sequence<pid> frontier;
    auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first);
    };
    auto make_pid = [&](node_id q) -> std::pair<node_id, double> {
      // Search for q in G.
      auto dist = distance(v[q]->coordinates.begin(), p_coords, d);
      return std::pair{q, dist};
    };
    int bits = std::ceil(std::log2(beamSize * beamSize));
    parlay::sequence<node_id> hash_table(1 << bits,
                                         std::numeric_limits<node_id>::max());

    // the frontier starts with the medoid
    frontier.push_back(make_pid(medoid->id));

    std::vector<pid> unvisited_frontier(beamSize);
    parlay::sequence<pid> new_frontier(2 * beamSize);
    unvisited_frontier[0] = frontier[0];
    bool not_done = true;

    // terminate beam search when the entire frontier has been visited
    while (not_done) {
      // the next node to visit is the unvisited frontier node that is closest
      // to p
      pid currentPid = unvisited_frontier[0];
      auto current_vtx = G.get_vertex(currentPid.first);

      auto g = [&](node_id a) {
        node_id loc = parlay::hash64_2(a) & ((1 << bits) - 1);
        if (hash_table[loc] == a) return false;
        hash_table[loc] = a;
        return true;
      };

      parlay::sequence<node_id> candidates;
      auto f = [&](node_id u, node_id v, empty_weight wgh) {
        if (g(v)) {
          candidates.push_back(v);
        }
        return true;
      };
      current_vtx.out_neighbors().foreach_cond(f);

      auto pairCandidates =
          parlay::map(candidates, [&](long c) { return make_pid(c); });
      auto sortedCandidates =
          parlay::unique(parlay::sort(pairCandidates, less));
      auto f_iter = std::set_union(
          frontier.begin(), frontier.end(), sortedCandidates.begin(),
          sortedCandidates.end(), new_frontier.begin(), less);
      size_t f_size = std::min<size_t>(beamSize, f_iter - new_frontier.begin());
      frontier =
          parlay::tabulate(f_size, [&](long i) { return new_frontier[i]; });
      visited.insert(
          std::upper_bound(visited.begin(), visited.end(), currentPid, less),
          currentPid);
      auto uf_iter =
          std::set_difference(frontier.begin(), frontier.end(), visited.begin(),
                              visited.end(), unvisited_frontier.begin(), less);
      not_done = uf_iter > unvisited_frontier.begin();
    }
#ifdef STATS
    total_visited.update_value(visited.size());
#endif
    return std::make_pair(frontier, parlay::to_sequence(visited));
  }

  // robustPrune routine as found in DiskANN paper.
  // The new candidate set is added to the supplied array (new_nghs).
  void robustPrune(tvec_point* p, node_id p_id,
                   parlay::sequence<pid> candidates, double alpha,
                   parlay::slice<node_id*, node_id*> new_nghs,
                   bool add = true) {
    auto less = [&](pid a, pid b) { return a.second < b.second; };

    // add out neighbors of p to the candidate set.
    if (add) {
      auto p_vertex = G.get_vertex(p_id);
      auto map_f = [&](node_id u, node_id neighbor_id, empty_weight weight) {
        candidates.push_back(std::make_pair(
            neighbor_id, distance(v[neighbor_id]->coordinates.begin(),
                                  p->coordinates.begin(), d)));
        return true;
      };
      p_vertex.out_neighbors().foreach_cond(map_f);
    }

    if (candidates.size() <= maxDeg) {
      uint32_t offset = 0;
      parlay::sort_inplace(candidates, less);
      for (size_t i = 0; i < candidates.size(); ++i) {
        node_id ngh = candidates[i].first;
        if (ngh != p_id && (i == 0 || ngh != candidates[i - 1].first)) {
          new_nghs[offset] = ngh;
          offset++;
        }
      }
      return;
    }

    // Sort the candidate set in reverse order according to distance from p.
    parlay::sort_inplace(candidates, less);

    size_t num_new = 0;

    size_t candidate_idx = 0;
    while (num_new <= maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      node_id p_star = candidates[candidate_idx].first;
      candidate_idx++;
      if (p_star == p->id || p_star == std::numeric_limits<node_id>::max()) {
        continue;
      }

      new_nghs[num_new] = p_star;
      num_new++;

      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        node_id p_prime = candidates[i].first;
        if (p_prime != std::numeric_limits<node_id>::max()) {
          float dist_starprime = distance(v[p_star]->coordinates.begin(),
                                          v[p_prime]->coordinates.begin(), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = std::numeric_limits<node_id>::max();
          }
        }
      }
    }
  }

  // wrapper to allow calling robustPrune on a sequence of candidates
  // that do not come with precomputed distances
  void robustPrune(Tvec_point<T>* p, node_id p_id,
                   parlay::sequence<node_id>& candidates, double alpha,
                   parlay::slice<node_id*, node_id*> new_nghs,
                   bool add = true) {
    parlay::sequence<pid> cc;
    auto p_vertex = G.get_vertex(p_id);
    node_id p_ngh_size = p_vertex.out_degree();

    cc.reserve(candidates.size() + p_ngh_size);
    for (size_t i = 0; i < candidates.size(); ++i) {
      cc.push_back(std::make_pair(
          candidates[i], distance(v[candidates[i]]->coordinates.begin(),
                                  p->coordinates.begin(), d)));
    }
    robustPrune(p, p_id, std::move(cc), alpha, new_nghs, add);
  }

  // special size function
  template <class G>
  static int size_of(parlay::slice<G*, G*> nbh) {
    int size = 0;
    size_t i = 0;
    while (nbh[i] != std::numeric_limits<node_id>::max() && i < nbh.size()) {
      size++;
      i++;
    }
    return size;
  }

  void batch_insert(parlay::sequence<node_id>& inserts, double base = 2,
                    double max_fraction = .02, bool random_order = true) {
    size_t n = v.size();
    size_t m = inserts.size();
    size_t max_batch_size = static_cast<size_t>(std::ceil(max_fraction * n));
    parlay::sequence<node_id> rperm;
    if (random_order)
      rperm = parlay::random_permutation<node_id>(static_cast<node_id>(m),
                                                  time(NULL));
    else
      rperm = parlay::tabulate(m, [&](node_id i) { return i; });
    auto shuffled_inserts =
        parlay::tabulate(m, [&](size_t i) { return inserts[rperm[i]]; });

    size_t floor = 0, ceiling = 0;
    uint32_t cnt_batch = 0;
    while (ceiling < m) {
      cnt_batch++;
      floor = ceiling;
      ceiling = std::min({m, size_t(floor * base) + 1, floor + max_batch_size});

      std::cout << "Start of batch insertion round: ceiling = " << ceiling
                << " floor = " << floor << std::endl;

      parlay::sequence<node_id> new_out = parlay::sequence<node_id>(
          maxDeg * (ceiling - floor), std::numeric_limits<node_id>::max());

      // search for each node starting from the medoid, then call
      // robustPrune with the visited list as its candidate set
      parlay::parallel_for(floor, ceiling, [&](size_t i) {
        size_t index = shuffled_inserts[i];

        parlay::sequence<pid> visited =
            (beam_search(v[index]->coordinates.begin(), beamSize)).second;
        auto new_nghs =
            parlay::make_slice(new_out.begin() + maxDeg * (i - floor),
                               new_out.begin() + maxDeg * (i + 1 - floor));
        robustPrune(v[index], index, visited, alpha, new_nghs);
      });
      // New neighbors of each point written into new_nbhs (above).
      // Not yet added to the graph G.

      // make each edge bidirectional by first adding each new edge
      // (i,j) to a sequence, then semisorting the sequence by key values
      auto to_flatten = parlay::tabulate(ceiling - floor, [&](size_t i) {
        node_id index = shuffled_inserts[i + floor];

        auto new_nghs = parlay::make_slice(new_out.begin() + maxDeg * i,
                                           new_out.begin() + maxDeg * (i + 1));

        auto edges = parlay::tabulate(size_of(new_nghs), [&](size_t j) {
          return std::make_pair((new_nghs)[j], index);
        });
        return edges;
      });

      auto KVs = parlay::tabulate(ceiling - floor, [&](size_t i) {
        node_id index = shuffled_inserts[i + floor];
        auto new_nghs = parlay::make_slice(new_out.begin() + maxDeg * i,
                                           new_out.begin() + maxDeg * (i + 1));
        size_t nghs_size = size_of(new_nghs);
        auto begin =
            (std::tuple<node_id, empty_weight>*)(new_out.begin() + maxDeg * i);
        auto tree = edge_tree(begin, begin + nghs_size);
        auto tree_ptr = tree.root;
        tree.root = nullptr;

        return std::make_tuple(index, tree_ptr);
      });

      G.insert_vertices_batch(KVs.size(), KVs.begin());
      std::cout << "After inserts, G.num_vertices() (max node_id) = "
                << G.num_vertices() << std::endl;

      // TODO: update the code below:
      auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
      auto reverse_KVs =
          parlay::sequence<std::tuple<node_id, edge_node*>>(grouped_by.size());

      std::cout << "Number of in-neighbors: " << grouped_by.size() << std::endl;
      // finally, add the bidirectional edges; if they do not make
      // the vertex exceed the degree bound, just add them to out_nbhs;
      // otherwise, use robustPrune on the vertex with user-specified alpha
      parlay::parallel_for(0, grouped_by.size(), [&](size_t j) {
        // std::cout << "here1" << std::endl;
        auto[index, candidates] = grouped_by[j];
        // TODO: simpler case when newsize <= maxDeg.
        parlay::sequence<node_id> new_out_2(
            maxDeg, std::numeric_limits<node_id>::max());
        auto output_slice =
            parlay::make_slice(new_out_2.begin(), new_out_2.begin() + maxDeg);
        // std::cout << "here2" << std::endl;
        robustPrune(v[index], index, candidates, alpha, output_slice);
        // std::cout << "here2" << std::endl;
        size_t deg = size_of(output_slice);
        // std::cout << "here3" << std::endl;
        auto begin = (std::tuple<node_id, empty_weight>*)new_out_2.begin();
        auto tree = edge_tree(begin, begin + deg);
        // std::cout << "here4" << std::endl;
        reverse_KVs[j] = {index, tree.root};
        tree.root = nullptr;
      });

      std::cout << "ReverseKVs.size = " << reverse_KVs.size() << std::endl;
      G.insert_vertices_batch(reverse_KVs.size(), reverse_KVs.begin());
    }
  }

  parlay::sequence<float> centroid_helper(slice_tvec a) {
    if (a.size() == 1) {
      parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
      for (unsigned i = 0; i < d; i++)
        centroid_coords[i] = static_cast<float>((a[0]->coordinates)[i]);
      return centroid_coords;
    } else {
      size_t n = a.size();
      parlay::sequence<float> c1;
      parlay::sequence<float> c2;
      parlay::par_do_if(n > 1000,
                        [&]() { c1 = centroid_helper(a.cut(0, n / 2)); },
                        [&]() { c2 = centroid_helper(a.cut(n / 2, n)); });
      parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
      for (unsigned i = 0; i < d; i++) {
        float result = (c1[i] + c2[i]) / 2;
        centroid_coords[i] = result;
      }
      return centroid_coords;
    }
  }

  // TODO(laxmand): Why fvec point here?
  tvec_point* medoid_helper(fvec_point* centroid, slice_tvec a) {
    if (a.size() == 1) {
      return a[0];
    } else {
      size_t n = a.size();
      tvec_point* a1;
      tvec_point* a2;
      parlay::par_do_if(
          n > 1000, [&]() { a1 = medoid_helper(centroid, a.cut(0, n / 2)); },
          [&]() { a2 = medoid_helper(centroid, a.cut(n / 2, n)); });
      float d1 =
          distance(centroid->coordinates.begin(), a1->coordinates.begin(), d);
      float d2 =
          distance(centroid->coordinates.begin(), a2->coordinates.begin(), d);
      if (d1 < d2)
        return a1;
      else
        return a2;
    }
  }

  // computes the centroid and then assigns the approx medoid as the point in v
  // closest to the centroid
  node_id find_approx_medoid() {
    parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v));
    fvec_point centroidp = Tvec_point<float>();
    centroidp.coordinates = parlay::make_slice(centroid);
    medoid = medoid_helper(&centroidp, parlay::make_slice(v));
    std::cout << "Medoid ID: " << medoid->id << std::endl;
    return medoid->id;
  }

  parlay::sequence<node_id> all_distances_sorted(node_id sample_id) {
    size_t n = G.num_vertices();
    auto sample_coord = v[sample_id]->coordinates.begin();
    auto dist_and_id =
        parlay::tabulate(n, [&](size_t i) -> std::pair<double, node_id> {
          auto i_coord = v[i]->coordinates.begin();
          return {distance(i_coord, sample_coord, d), i};
        });
    parlay::sort_inplace(dist_and_id);
    return parlay::map(dist_and_id,
                       [&](const auto& pair) { return pair.second; });
  }

  void compute_self_recall(size_t num_samples = 1000) {
    size_t n = G.num_vertices();
    auto ids = parlay::tabulate(num_samples, [&](node_id i) -> node_id {
      return parlay::hash32(i) % n;
    });
    std::vector<int> ks = {1, 2, 4, 8, 16, 32, 64, 128, 256};
    auto counts = parlay::sequence<std::vector<int>>::from_function(
        num_samples, [&](size_t i) { return std::vector<int>(ks.size()); });
    parlay::parallel_for(0, num_samples, [&](size_t i) {
      auto sample_id = ids[i];
      auto ground_truth_nn = all_distances_sorted(sample_id);
      for (size_t j = 0; j < ks.size(); ++j) {
        size_t k = ks[j];
        auto neighbors = query(v[sample_id]->coordinates.begin(), k,
                               std::max((size_t)(2 * k + 1), (size_t)maxDeg));
        auto tab =
            std::set(ground_truth_nn.begin(), ground_truth_nn.begin() + k);
        size_t hits = 0;
        for (auto ngh : neighbors) {
          if (tab.find(ngh) != tab.end()) ++hits;
        }
        counts[i][j] = hits;
      }
    });
    for (size_t j = 0; j < ks.size(); ++j) {
      size_t k = ks[j];
      auto ds = parlay::delayed_seq<size_t>(
          num_samples, [&](size_t i) { return counts[i][j]; });
      auto tot = parlay::reduce(ds);
      std::cout << k << "@" << k << ": "
                << (static_cast<double>(tot) / (k * num_samples)) << std::endl;
    }
  }
};
