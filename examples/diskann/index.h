#pragma once

#include <cpam/cpam.h>
#include <pam/pam.h>
// using namespace cpam;
#include <pam/get_time.h>
#include <pam/parse_command_line.h>

#include "types.h"
#include "util/NSGDist.h"
#include "util/distance.h"

// Want to store a *directed* graph (without in-edges) using CPAM
// format. No augmentation needed.
// Vertices stored in a CPAM map, with modest block size. Each vertex
// object just points to its incident edges.

// Graph stores adjacency information (neighbors).

template <typename T>
struct ann_node_data {
  parlay::sequence<uint32_t> neighbor_ids;
  using tvec_point = Tvec_point<T>;

  tvec_point* our_point;

  // Other methods on a node? Do we store distances?
};

template <typename T>
struct ann_node_entry {
  using key_t = uint32_t;  // integer id
  using val_t = ann_node_data<T>;
  static inline bool comp(key_t a, key_t b) { return a < b; }
  using entry_t = std::tuple<key_t, val_t>;
};

template <typename T>
struct ann_graph {
  using ann_node_tree = cpam::pam_map<ann_node_entry<T>>;
  // Useful definitions
  using ann_node = typename ann_node_tree::node;
  using ann_node_gc = typename ann_node_tree::GC;
};

template <typename T>
struct knn_index {
  size_t maxDeg;
  size_t beamSize;
  // TODO(laxmand, magdalen): rename?
  double r2_alpha;  // alpha parameter for round 2 of robustPrune
  size_t d;

  using tvec_point = Tvec_point<T>;
  using fvec_point = Tvec_point<float>;

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
    // find the medoid, which each beamSearch will begin from
    find_approx_medoid(v);
    batch_insert(inserts, v, 2, .02, false);
    // build_index_inner(v, r2_alpha, 2, .02);
  }

 private:

  template <typename T>
  std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
      Tvec_point<T>* p, parlay::sequence<Tvec_point<T>*>& v,
      Tvec_point<T>* medoid, int beamSize, unsigned d) {
    // initialize data structures
    std::vector<pid> visited;
    parlay::sequence<pid> frontier;
    auto less = [&](pid a, pid b) {
        return a.second < b.second || (a.second == b.second && a.first < b.first); };
    auto make_pid = [&] (int q) {
      return std::pair{q, distance(v[q]->coordinates.begin(), p->coordinates.begin(), d)};};
    int bits = std::ceil(std::log2(beamSize*beamSize));
    parlay::sequence<int> hash_table(1 << bits, -1);

    // the frontier starts with the medoid
    frontier.push_back(make_pid(medoid->id));

    std::vector<pid> unvisited_frontier(beamSize);
    parlay::sequence<pid> new_frontier(2*beamSize);
    unvisited_frontier[0] = frontier[0];
    bool not_done = true;

    // terminate beam search when the entire frontier has been visited
    while (not_done) {
      // the next node to visit is the unvisited frontier node that is closest to p
      pid currentPid = unvisited_frontier[0];
      Tvec_point<T>* current = v[currentPid.first];
      auto g = [&] (int a) {
  	       int loc = parlay::hash64_2(a) & ((1 << bits) - 1);
  	       if (hash_table[loc] == a) return false;
  	       hash_table[loc] = a;
  	       return true;};

      // TODO: Neighborhood access:
      auto candidates = parlay::filter(current->out_nbh.cut(0,size_of(current->out_nbh)), g);

      auto pairCandidates = parlay::map(candidates, [&] (long c) {return make_pid(c);});
      auto sortedCandidates = parlay::unique(parlay::sort(pairCandidates, less));
      auto f_iter = std::set_union(frontier.begin(), frontier.end(),
  				 sortedCandidates.begin(), sortedCandidates.end(),
  				 new_frontier.begin(), less);
      size_t f_size = std::min<size_t>(beamSize, f_iter - new_frontier.begin());
      frontier = parlay::tabulate(f_size, [&] (long i) {return new_frontier[i];});
      visited.insert(std::upper_bound(visited.begin(), visited.end(), currentPid, less), currentPid);
      auto uf_iter = std::set_difference(frontier.begin(), frontier.end(),
  				 visited.begin(), visited.end(),
  				 unvisited_frontier.begin(), less);
      not_done = uf_iter > unvisited_frontier.begin();
    }
    return std::make_pair(frontier, parlay::to_sequence(visited));
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

      std::cout << "here3: ceiling = " << ceiling << " floor = " << floor << std::endl;

      parlay::sequence<int> new_out =
          parlay::sequence<int>(maxDeg * (ceiling - floor), -1);

//      // search for each node starting from the medoid, then call
//      // robustPrune with the visited list as its candidate set
//      parlay::parallel_for(floor, ceiling, [&](size_t i) {
//        size_t index = shuffled_inserts[i];
//        v[index]->new_nbh =
//            parlay::make_slice(new_out.begin() + maxDeg * (i - floor),
//                               new_out.begin() + maxDeg * (i + 1 - floor));
//        parlay::sequence<pid> visited =
//            (beam_search(v[index], v, medoid, beamSize, d)).second;
//        if (report_stats) v[index]->cnt = visited.size();
//        robustPrune(v[index], visited, v, r2_alpha);
//      });
//
//      // make each edge bidirectional by first adding each new edge
//      //(i,j) to a sequence, then semisorting the sequence by key values
//      // std::cout << "here3" << std::endl;
//      auto to_flatten = parlay::tabulate(ceiling - floor, [&](size_t i) {
//        int index = shuffled_inserts[i + floor];
//        auto edges =
//            parlay::tabulate(size_of(v[index]->new_nbh), [&](size_t j) {
//              return std::make_pair((v[index]->new_nbh)[j], index);
//            });
//        return edges;
//      });
//      parlay::parallel_for(floor, ceiling, [&](size_t i) {
//        synchronize(v[shuffled_inserts[i]]);
//      });
//      auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
//      // finally, add the bidirectional edges; if they do not make
//      // the vertex exceed the degree bound, just add them to out_nbhs;
//      // otherwise, use robustPrune on the vertex with user-specified alpha
//      // std::cout << "here4" << std::endl;
//      parlay::parallel_for(0, grouped_by.size(), [&](size_t j) {
//        auto[index, candidates] = grouped_by[j];
//        int newsize = candidates.size() + size_of(v[index]->out_nbh);
//        if (newsize <= maxDeg) {
//          for (const int& k : candidates) add_nbh(k, v[index]);
//        } else {
//          parlay::sequence<int> new_out_2(maxDeg, -1);
//          v[index]->new_nbh =
//              parlay::make_slice(new_out_2.begin(), new_out_2.begin() + maxDeg);
//          robustPrune(v[index], candidates, v, r2_alpha);
//          synchronize(v[index]);
//        }
//      });
//      // std::cout << "here5" << std::endl;
//      inc += 1;

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
  void find_approx_medoid(parlay::sequence<Tvec_point<T>*>& v) {
    parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v));
    fvec_point centroidp = Tvec_point<float>();
    centroidp.coordinates = parlay::make_slice(centroid);
    medoid = medoid_helper(&centroidp, parlay::make_slice(v));
    std::cout << "Medoid ID: " << medoid->id << std::endl;
  }
};
