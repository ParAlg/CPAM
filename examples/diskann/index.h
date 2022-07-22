#pragma once

#include <cpam/cpam.h>
#include <pam/pam.h>
// using namespace cpam;
#include <pam/get_time.h>
#include <pam/parse_command_line.h>

#include "types.h"
#include "distance.h"

// Want to store a *directed* graph (without in-edges) using CPAM
// format. No augmentation needed.
// Vertices stored in a CPAM map, with modest block size. Each vertex
// object just points to its incident edges.

// Graph stores adjacency information (neighbors).

template<typename T>
void clear(Tvec_point<T>* p) {
	for(int j=0; j<p->out_nbh.size(); j++) p->out_nbh[j] = -1;
}

template<typename T>
void clear(parlay::sequence<Tvec_point<T>*> &v) {
	size_t n = v.size();
	parlay::parallel_for(0, n, [&] (size_t i){
		for(int j=0; j<v[i]->out_nbh.size(); j++) v[i]->out_nbh[j] = -1;
	});
}

template <typename T>
struct ann_node_data {
  parlay::sequence<uint32_t> neighbor_ids;
  using tvec_point = Tvec_point<T>;

  // TODO: pointer or index?
  tvec_point* our_point;

  // Other methods on a node? Do we store distances?
};

template <typename T>
struct ann_node_entry {
  using key_t = uint32_t; // integer id
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

  // TODO(laxmand): why use int for inserts?
	void build_index(parlay::sequence<Tvec_point<T>*> &v, parlay::sequence<int> inserts) {
		clear(v);
		//find the medoid, which each beamSearch will begin from
		find_approx_medoid(v);
		// batch_insert(inserts, v, 2, .02, false);
		// build_index_inner(v, r2_alpha, 2, .02);
	}

 private:
  parlay::sequence<float> centroid_helper(slice_tvec a) {
    if (a.size() == 1) {
      parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
      for (int i = 0; i < d; i++)
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
      for (int i = 0; i < d; i++) {
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
    // TODO(laxmand): ?? Fix Distance?
      float d1, d2;
//      float d1 =
//          distance(centroid->coordinates.begin(), a1->coordinates.begin(), d);
//      float d2 =
//          distance(centroid->coordinates.begin(), a2->coordinates.begin(), d);
      if (d1 < d2)
        return a1;
      else
        return a2;
    }
  }

  // computes the centroid and then assigns the approx medoid as the point in v
  // closest to the centroid
  void find_approx_medoid(parlay::sequence<Tvec_point<T>*>& v) {
    size_t n = v.size();
    parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v));
    fvec_point centroidp = Tvec_point<float>();
    centroidp.coordinates = parlay::make_slice(centroid);
    medoid = medoid_helper(&centroidp, parlay::make_slice(v));
    std::cout << "Medoid ID: " << medoid->id << std::endl;
  }

};
