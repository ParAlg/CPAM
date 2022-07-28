#pragma once

#include <cpam/cpam.h>
#include <pam/pam.h>
// using namespace cpam;
#include <pam/get_time.h>
#include <pam/parse_command_line.h>

#include "../graphs/aspen/aspen.h"
#include "types.h"
#include "util/NSGDist.h"
#include "util/distance.h"

bool stats = false;

template <typename T>
struct knn_index {
  // Pointers to the coordinates for each point.
  parlay::sequence<Tvec_point<T>*>& v;

  size_t maxDeg;
  size_t beamSize;
  // TODO: rename? are we using round 2?
  double alpha;  // alpha parameter for round 2 of robustPrune
  size_t d;

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

  using pid = std::pair<int, float>;
  using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
  using index_pair = std::pair<int, int>;
  using slice_idx = decltype(make_slice(parlay::sequence<index_pair>()));

  knn_index(parlay::sequence<Tvec_point<T>*>& v, size_t maxDeg, size_t beamSize,
            double Alpha, size_t dim)
      : v(v), maxDeg(maxDeg), beamSize(beamSize), alpha(Alpha), d(dim) {
    std::cout << "Initialized knn_index with maxDeg = " << maxDeg
              << " beamSize = " << beamSize << " alpha = " << alpha
              << " dim = " << dim << std::endl;
  }

  void build_index(parlay::sequence<int> inserts) {
    // Find the medoid, which each beamSearch will begin from.
    node_id medoid_id = find_approx_medoid();
    std::cout << "medoid coordinates = " << medoid->coordinates.begin()
              << std::endl;
    G.insert_vertex_inplace(medoid_id, nullptr);
    batch_insert(inserts, 2, .02, false);
    std::cout << "G.num_vertices = " << G.num_vertices()
              << " num_edges = " << G.num_edges() << std::endl;

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
    for (int p : deletes) {
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

  // void consolidate_deletes(){

  // }

 private:
  std::set<int> delete_set;
  // p_coords: query vector coordinates
  // v: database of vectors
  // medoid: "root" of the proximity graph
  // beamSize: (similar to ef)
  // d: dimensionality of the indexed vectors
  std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
      T* p_coords, int beamSize) {
    // initialize data structures
    // pid = std::pair<int, float>
    std::vector<pid> visited;
    parlay::sequence<pid> frontier;
    auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first);
    };
    auto make_pid = [&](int q) -> std::pair<int, double> {
      // Search for q in G.
      auto dist = distance(v[q]->coordinates.begin(), p_coords, d);
      return std::pair{q, dist};
    };
    int bits = std::ceil(std::log2(beamSize * beamSize));
    parlay::sequence<int> hash_table(1 << bits, -1);

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

      auto g = [&](int a) {
        int loc = parlay::hash64_2(a) & ((1 << bits) - 1);
        if (hash_table[loc] == a) return false;
        hash_table[loc] = a;
        return true;
      };

      // TODO: reserve?
      parlay::sequence<int> candidates;
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
    return std::make_pair(frontier, parlay::to_sequence(visited));
  }

  // robustPrune routine as found in DiskANN paper.
  // The new candidate set is added to the supplied array (new_nghs).
  void robustPrune(tvec_point* p, node_id p_id,
                   parlay::sequence<pid> candidates, double alpha,
                   parlay::slice<int*, int*> new_nghs, bool add = true) {
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

    // TODO: Is this snippet correct? Do we care bout duplicates?
    if (candidates.size() <= maxDeg) {
      for (size_t i = 0; i < candidates.size(); ++i) {
        new_nghs[i] = candidates[i].first;
      }
      return;
    }

    // Sort the candidate set in reverse order according to distance from p.
    auto less = [&](pid a, pid b) { return a.second < b.second; };
    parlay::sort_inplace(candidates, less);

    parlay::sequence<int> new_nbhs = parlay::sequence<int>();

    size_t candidate_idx = 0;
    while (new_nbhs.size() <= maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      int p_star = candidates[candidate_idx].first;
      candidate_idx++;
      if (p_star == p->id || p_star == -1) {
        continue;
      }

      new_nbhs.push_back(p_star);

      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        int p_prime = candidates[i].first;
        if (p_prime != -1) {
          float dist_starprime = distance(v[p_star]->coordinates.begin(),
                                          v[p_prime]->coordinates.begin(), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = -1;
          }
        }
      }
    }
    for (size_t i = 0; i < new_nbhs.size(); ++i) {
      new_nghs[i] = new_nbhs[i];  // change names
    }
  }

  // wrapper to allow calling robustPrune on a sequence of candidates
  // that do not come with precomputed distances
  void robustPrune(Tvec_point<T>* p, node_id p_id,
                   parlay::sequence<int>& candidates, double alpha,
                   parlay::slice<int*, int*> new_nghs, bool add = true) {
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
    while (nbh[i] != -1 && i < nbh.size()) {
      size++;
      i++;
    }
    return size;
  }

  void batch_insert(parlay::sequence<int>& inserts, double base = 2,
                    double max_fraction = .02, bool random_order = true) {
    std::cout << "here" << std::endl;

    size_t n = v.size();
    size_t m = inserts.size();
    size_t inc = 0;
    size_t count = 0;
    size_t max_batch_size = static_cast<size_t>(max_fraction * n);
    parlay::sequence<int> rperm;
    if (random_order)
      rperm = parlay::random_permutation<int>(static_cast<int>(m), time(NULL));
    else
      rperm = parlay::tabulate(m, [&](int i) { return i; });
    auto shuffled_inserts =
        parlay::tabulate(m, [&](size_t i) { return inserts[rperm[i]]; });
    std::cout << "here1" << std::endl;

    size_t floor;
    size_t ceiling;
    while (count < m) {
      if (pow(base, inc) <= max_batch_size) {
        floor = static_cast<size_t>(pow(base, inc)) - 1;
        ceiling = std::min(static_cast<size_t>(pow(base, inc + 1)), m) - 1;
        count = std::min(static_cast<size_t>(pow(base, inc + 1)), m) - 1;
      } else {
        floor = count;
        ceiling = std::min(count + static_cast<size_t>(max_batch_size), m);
        count += static_cast<size_t>(max_batch_size);
      }

      std::cout << "Start of batch insertion round: ceiling = " << ceiling
                << " floor = " << floor << std::endl;

      parlay::sequence<int> new_out =
          parlay::sequence<int>(maxDeg * (ceiling - floor), -1);

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
      // std::cout << "here3" << std::endl;
      auto to_flatten = parlay::tabulate(ceiling - floor, [&](size_t i) {
        int index = shuffled_inserts[i + floor];

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

      std::cout << "Number of vertex insertions for batch: " << inc << " : "
                << KVs.size() << std::endl;
      G.insert_vertices_batch(KVs.size(), KVs.begin());
      std::cout << "After inserts, G.num_vertices() (max node_id) = "
                << G.num_vertices() << std::endl;

      // TODO: update the code below:
      auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
      auto reverse_KVs =
          parlay::sequence<std::tuple<node_id, edge_node*>>(grouped_by.size());

      // finally, add the bidirectional edges; if they do not make
      // the vertex exceed the degree bound, just add them to out_nbhs;
      // otherwise, use robustPrune on the vertex with user-specified alpha
      parlay::parallel_for(0, grouped_by.size(), [&](size_t j) {
        auto[index, candidates] = grouped_by[j];
        // TODO: simpler case when newsize <= maxDeg.
        parlay::sequence<int> new_out_2(maxDeg, -1);
        auto output_slice =
            parlay::make_slice(new_out_2.begin(), new_out_2.begin() + maxDeg);
        robustPrune(v[index], index, candidates, alpha, output_slice);
        size_t deg = size_of(output_slice);
        auto begin = (std::tuple<node_id, empty_weight>*)new_out_2.begin();
        auto tree = edge_tree(begin, begin + deg);
        reverse_KVs[j] = {index, tree.root};
        tree.root = nullptr;
      });

      std::cout << "ReverseKVs.size = " << reverse_KVs.size() << std::endl;
      G.insert_vertices_batch(reverse_KVs.size(), reverse_KVs.begin());

      inc += 1;
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
