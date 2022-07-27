#pragma once

#include <cpam/cpam.h>
#include <pam/get_time.h>
#include <pam/pam.h>
#include <pam/parse_command_line.h>

//#include "types.h"

//// Set
//struct edge_entry {
//  using key_t = node_id;
//  static inline bool comp(key_t a, key_t b) { return a < b; }
//};
//
//using edge_tree = cpam::pam_set<edge_entry>;
//using edge_node = typename edge_tree::node;
//
//// Map
//template <typename T>
//struct node_entry {
//  using key_t = node_id;           // Integer id of this node.
//  using val_t = edge_node*;    // Set storing the neighbors.
//  static inline bool comp(key_t a, key_t b) { return a < b; }
//  using entry_t = std::tuple<key_t, val_t>;
//};

//// T is the type of the underlying vector data.
//// B is the block size used in the CPAM map.
//template <typename T, size_t B = 128>
//struct ann_graph {
//  using node_tree = cpam::pam_map<node_entry<T>, B>;
//  using ET = typename node_entry<T>::entry_t;
//
//  node_tree V;
//
//  ann_graph() {}
//
//  // Note that vertices are not stored but are just wrappers around
//  // the underlying node_entry with some helper functions.
//  struct readonly_node {
//    node_id id;
//    edge_node* node;
//
//    // Return the out-degree of the node in the directed graph.
//    size_t out_degree() const { return edge_tree::size(node); }
//
//    template <class G>
//    parlay::sequence<node_id> filter(const G& g) const {
//      // TODO
////      auto s = parlay::make_slice(node.neighbors, node.neighbors + node.num_neighbors);
////      return parlay::filter(s, g);
//    }
//
//    template <class F>
//    void map(const F& f) const {
//      auto map_f = [&](const auto& et, size_t i) {
//        f(et);
//      };
//      edge_tree tree;
//      tree.root = node;
//      tree.foreach_index(tree, map_f);
//      tree.root = nullptr;
//    }
//
//    read_only_node() : id(std::numeric_limits<node_id>::max()), node(nullptr) {}
//    readonly_node(node_id id, edge_node* node) : id(id), node(node) {}
//  };
//
//
//  // For now, going to use G as an in-place data structure. We'll add
//  // the versioning DS in the next step.
//  void insert_empty_node(node_id id) {
//    auto et = std::make_tuple(id, nullptr);
//    V.insert(et);
//  }
//
//  void batch_insert_inplace(parlay::sequence<ET> entries) {
//    auto replace = [&] (edge_node* prev_node, edge_node* new_node) {
//      assert(false);  // TODO: In general this case can occur---need to dealloc prev_node.
//      return new_node;
//    };
//    V = node_tree::multi_insert_sorted(std::move(V), parlay::make_slice(entries), replace);
//  }
//
//  vertex get_node(node_id id) {
//    auto opt = V.find(id);
//    if (opt.has_value()) {
//      const auto& in_opt = *opt;
//      // auto ref_cnt = edge_tree::Tree::ref_cnt(in_opt);
//      // assert(ref_cnt == 1);
//      return readonly_node(id, in_opt);
//    }
//    return readonly_node(id, nullptr);
//  }
//
//};
