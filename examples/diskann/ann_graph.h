#pragma once

#include <cpam/cpam.h>
#include <pam/get_time.h>
#include <pam/pam.h>
#include <pam/parse_command_line.h>

#include "types.h"

template <typename T>
struct ann_node_data {
  // The start of the coordinates for this point. The dimensionality d
  // is known by the index, so there's no need to redundantly store it
  // per-point.
  T* coordinates;

  // Ids of the neighboring nodes of this point in the underlying
  // proximity graph.
  node_id num_neighbors;
  node_id* neighbors;

  ann_node_data()
      : coordinates(nullptr), num_neighbors(0), neighbors(nullptr) {}
  ann_node_data(T* coordinates, node_id num_neighbors, node_id* neighbors)
      : coordinates(coordinates),
        num_neighbors(num_neighbors),
        neighbors(neighbors) {}
};

template <typename T>
struct ann_node_entry {
  using key_t = node_id;           // Integer id of this node.
  using val_t = ann_node_data<T>;  // Coordinates and neighbors of this node.
  static inline bool comp(key_t a, key_t b) { return a < b; }
  using entry_t = std::tuple<key_t, val_t>;
};

// T is the type of the underlying vector data.
// B is the block size used in the CPAM map.
template <typename T, size_t B = 128>
struct ann_graph {
  using ann_node_tree = cpam::pam_map<ann_node_entry<T>, B>;
  using ann_node = ann_node_data<T>;
  using ann_entry = typename ann_node_entry<T>::entry_t;

  ann_node_tree V;

  ann_graph() {}

  // Note that vertices are not stored but are just wrappers around
  // the underlying ann_node_entry with some helper functions.
  struct readonly_node {
    size_t id;
    const ann_node& node;

    // Return the out-degree of the node in the directed graph.
    // Note that in-degree is not supported.
    size_t out_degree() const { return node.neighbor_ids.size(); }

    template <class G>
    parlay::sequence<node_id> filter(const G& g) const {
      auto s = parlay::make_slice(node.neighbors, node.neighbors + node.num_neighbors);
      return parlay::filter(s, g);
    }

    template <class F>
    void map(const F& f) const {
      auto s = parlay::make_slice(node.neighbors, node.neighbors + node.num_neighbors);
      for (const auto& x : s) {
        f(x);
      }
    }

    readonly_node(node_id id, const ann_node& node) : id(id), node(node) {}
  };

  // For now, going to use G as an in-place data structure. We'll add
  // the versioning DS in the next step.
  void insert_node_inplace(node_id id, const ann_node& data) {
    auto et = std::make_tuple(id, data);
    V.insert(et);
  }

  void batch_insert_inplace(parlay::sequence<ann_entry> entries) {
    auto replace = [&] (ann_node prev_node, ann_node new_node) {
      assert(false);  // TODO: In general this case can occur---need to dealloc prev_node.
      return new_node;
    };
    V = multi_insert_sorted(std::move(V), entries, replace);
  }

  std::optional<readonly_node> get_node(node_id id) {
    auto opt = V.find(id);
    if (opt.has_value()) {
      const ann_node& in_opt = *opt;
      // auto ref_cnt = edge_tree::Tree::ref_cnt(in_opt);
      // assert(ref_cnt == 1);
      return readonly_node(id, in_opt);
    }
    return {};
  }

  T* get_coordinates(node_id id) {
    auto opt = V.find(id);
    if (opt.has_value()) {
      return opt->coordinates;
    }
    return nullptr;
  }

};
