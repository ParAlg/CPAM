#pragma once

#include <cpam/cpam.h>
#include <pam/pam.h>
// using namespace cpam;
#include <pam/get_time.h>
#include <pam/parse_command_line.h>

#include "ann_graph.h"
#include "types.h"
#include "util/NSGDist.h"
#include "util/distance.h"

template <typename T>
struct knn_index {
  size_t maxDeg;
  size_t beamSize;
  // TODO(laxmand, magdalen): rename?
  double r2_alpha;  // alpha parameter for round 2 of robustPrune
  size_t d;

  using ANNGraph = ann_graph<T>;
  ANNGraph G;

  using tvec_point = Tvec_point<T>;
  using fvec_point = Tvec_point<float>;
  using readonly_node = typename ANNGraph::readonly_node;

  tvec_point* medoid;

  using pid = std::pair<int, float>;
  using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));
  using index_pair = std::pair<int, int>;
  using slice_idx = decltype(make_slice(parlay::sequence<index_pair>()));
  using fine_sequence = parlay::sequence<int>;

  knn_index(size_t maxDeg, size_t beamSize, double alpha, size_t dim)
      : maxDeg(maxDeg), beamSize(beamSize), r2_alpha(alpha), d(dim) {}

  void build_index(parlay::sequence<Tvec_point<T>*>& v,
                   parlay::sequence<int> inserts) {
    // Find the medoid, which each beamSearch will begin from.
    node_id medoid_id = find_approx_medoid(v);
    auto node_data = ann_node_data<T>(medoid->coordinates.begin(), 0, nullptr);
    G.insert_node_inplace(medoid_id, node_data);
    auto x = G.get_node(medoid_id);
    std::cout << "here: " << x.has_value() << std::endl;
    // batch_insert(inserts, v, 2, .02, false);
  }

 private:
  // p_coords: query vector coordinates
  // v: database of vectors
  // medoid: "root" of the proximity graph
  // beamSize: (similar to ef)
  // d: dimensionality of the indexed vectors
  std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
      T* p_coords, Tvec_point<T>* medoid, int beamSize, unsigned d) {
    // initialize data structures
    // pid = std::pair<int, float>
    std::vector<pid> visited;
    parlay::sequence<pid> frontier;
    auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first);
    };
    auto make_pid = [&](int q) {
      // Search for q in G.
      auto v_q_opt = G.get_node(q);
      assert(v_q.has_value());
      auto v_q = *v_q_opt;
      return std::pair{q, distance(v_q.coordinates, p_coords, d)};
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
      auto current_opt = G.get_node(currentPid.first);
      assert(current_opt.has_value());
      const auto& current = *current;

      auto g = [&](int a) {
        int loc = parlay::hash64_2(a) & ((1 << bits) - 1);
        if (hash_table[loc] == a) return false;
        hash_table[loc] = a;
        return true;
      };

      // TODO: Neighborhood access:
      auto candidates = current.filter(g);

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
                   parlay::sequence<pid> candidates,
                   parlay::sequence<tvec_point*>& v, double alpha,
                   parlay::slice<int*, int*> new_nghs, bool add = true) {
    // add out neighbors of p to the candidate set.
    if (add) {
      auto p_node_opt = G.get_node(p_id);
      if (p_node_opt.has_value()) {
        auto map_f = [&](node_id neighbor_id) {
          candidates.push_back(std::make_pair(
              neighbor_id,
              distance(G.get_coordinates(neighbor_id), p->coordinates.begin(), d)));
        };
        p_node_opt->map(map_f);
      }
    }

    // Sort the candidate set in reverse order according to distance from p.
    auto less = [&](pid a, pid b) { return a.second < b.second; };
    parlay::sort_inplace(candidates, less);

    fine_sequence new_nbhs = fine_sequence();

    size_t candidate_idx = 0;
    while (new_nbhs.size() <= maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      int p_star = candidates[candidate_idx].first;
      candidate_idx++;
      if (p_star == p->id || p_star == -1) {
        continue;
      }

      new_nbhs.push_back(p_star);

      // TODO: easy optimization: do all get_coordinate calls up-front
      // and cache them.
      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        int p_prime = candidates[i].first;
        if (p_prime != -1) {
          float dist_starprime = distance(G.get_coordinates(p_star),
                                          G.get_coordinates(p_prime), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = -1;
          }
        }
      }
    }
    for (size_t i=0; i<new_nbhs.size(); ++i) {
      new_nghs[i] = new_nbhs[i];  // change names
    }
  }

  // wrapper to allow calling robustPrune on a sequence of candidates
  // that do not come with precomputed distances
  void robustPrune(Tvec_point<T>* p, node_id p_id,
                   parlay::sequence<int> candidates,
                   parlay::sequence<Tvec_point<T>*>& v, double alpha,
                   parlay::slice<int*, int*> new_nghs, bool add = true) {
    parlay::sequence<pid> cc;
    auto p_node_opt = G.get_node(p_id);
    node_id p_ngh_size = p_node_opt.has_value ? *p_node_opt.out_degree() : 0;

    cc.reserve(candidates.size() + p_ngh_size);
    for (size_t i = 0; i < candidates.size(); ++i) {
      auto v_q_opt = G.get_node(candidates[i]);
      assert(v_q.has_value());
      auto v_q = *v_q_opt;
      cc.push_back(std::make_pair(
          candidates[i], distance(v_q.coordinates, p->coordinates.begin(), d)));
    }
    robustPrune(p, p_id, std::move(cc), v, alpha, new_nghs, add);
  }

  //special size function
  static int size_of(parlay::slice<T*, T*> nbh) const {
  	int size = 0;
  	int i=0;
  	while(nbh[i] != -1 && i<nbh.size()) {size++; i++;}
  	return size;
  }

  void batch_insert(parlay::sequence<int>& inserts,
                    parlay::sequence<Tvec_point<T>*>& v, double base = 2,
                    double max_fraction = .02, bool random_order = true) {
    for (int p : inserts) {
      if (p < 0 || p > (int)v.size()) {
        std::cout << "ERROR: invalid or already inserted point " << p
                  << " given to batch_insert" << std::endl;
        abort();
      }
    }
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

    while (count < m) {
      size_t floor;
      size_t ceiling;
      if (pow(base, inc) <= max_batch_size) {
        floor = static_cast<size_t>(pow(base, inc)) - 1;
        ceiling = std::min(static_cast<size_t>(pow(base, inc + 1)), m) - 1;
        count = std::min(static_cast<size_t>(pow(base, inc + 1)), m) - 1;
        // std::cout << "here2" << std::endl;
      } else {
        floor = count;
        ceiling = std::min(count + static_cast<size_t>(max_batch_size), m) - 1;
        count += static_cast<size_t>(max_batch_size);
      }

      std::cout << "here3: ceiling = " << ceiling << " floor = " << floor
                << std::endl;

      parlay::sequence<int> new_out =
          parlay::sequence<int>(maxDeg * (ceiling - floor), -1);

      // search for each node starting from the medoid, then call
      // robustPrune with the visited list as its candidate set
      parlay::parallel_for(floor, ceiling, [&](size_t i) {
        size_t index = shuffled_inserts[i];
        // v[index]->new_nbh =
        //    parlay::make_slice(new_out.begin() + maxDeg * (i - floor),
        //                       new_out.begin() + maxDeg * (i + 1 - floor));

        parlay::sequence<pid> visited =
            (beam_search(v[index]->coordinates.begin(), medoid, beamSize, d))
                .second;
        // if (report_stats) v[index]->cnt = visited.size();
        auto new_nghs =
            parlay::make_slice(new_out.begin() + maxDeg * (i - floor),
                               new_out.begin() + maxDeg * (i + 1 - floor));
        robustPrune(v[index], index, visited, v, r2_alpha, new_nghs);
      });
      // New neighbors of each point written into new_nbhs (above).
      // Not yet added to the graph G.

      // make each edge bidirectional by first adding each new edge
      //(i,j) to a sequence, then semisorting the sequence by key values
      // std::cout << "here3" << std::endl;
      auto to_flatten = parlay::tabulate(ceiling - floor, [&](size_t i) {
        int index = shuffled_inserts[i + floor];
        auto edges =
            parlay::tabulate(size_of(v[index]->new_nbh), [&](size_t j) {
              return std::make_pair((v[index]->new_nbh)[j], index);
            });
        return edges;
      });
      parlay::parallel_for(floor, ceiling, [&](size_t i) {
        synchronize(v[shuffled_inserts[i]]);
      });
      auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
      // finally, add the bidirectional edges; if they do not make
      // the vertex exceed the degree bound, just add them to out_nbhs;
      // otherwise, use robustPrune on the vertex with user-specified alpha
      // std::cout << "here4" << std::endl;
      parlay::parallel_for(0, grouped_by.size(), [&](size_t j) {
        auto[index, candidates] = grouped_by[j];
        int newsize = candidates.size() + size_of(v[index]->out_nbh);
        if (newsize <= maxDeg) {
          for (const int& k : candidates) add_nbh(k, v[index]);
        } else {
          parlay::sequence<int> new_out_2(maxDeg, -1);
          v[index]->new_nbh =
              parlay::make_slice(new_out_2.begin(), new_out_2.begin() + maxDeg);
          robustPrune(v[index], candidates, v, r2_alpha);
          synchronize(v[index]);
        }
      });

      // std::cout << "here5" << std::endl;
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
  node_id find_approx_medoid(parlay::sequence<Tvec_point<T>*>& v) {
    parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v));
    fvec_point centroidp = Tvec_point<float>();
    centroidp.coordinates = parlay::make_slice(centroid);
    medoid = medoid_helper(&centroidp, parlay::make_slice(v));
    std::cout << "Medoid ID: " << medoid->id << std::endl;
    return medoid->id;
  }
};
