#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>
#include <random>

#include <parlay/io.h>
#include <parlay/primitives.h>
#include <parlay/random.h>


#include <pam/get_time.h>
#include <pam/parse_command_line.h>
#include "../../../graphs/aspen/aspen.h"
#include "../../index.h"

// using node_id = long int; 

struct test_graph{
  struct empty_weight {};
  using Graph = aspen::symmetric_graph<empty_weight>;
  using vertex = typename Graph::vertex;
  using edge_tree = typename Graph::edge_tree;
  using edge_node = typename Graph::edge_node;
  using vertex_tree = typename Graph::vertex_tree;
  using vertex_node = typename Graph::vertex_node;
  aspen::versioned_graph<Graph> VG;

  void insert(size_t n, parlay::sequence<node_id> inserts, parlay::sequence<node_id> prev_inserts){
    auto S = VG.acquire_version();
    std::cout << "Acquired version with timestamp " << S.timestamp
              << " address = " << ((size_t)S.get_root())
              << " ref_cnt = " << S.get_ref_cnt() << std::endl;

    parlay::sequence<node_id> init_inserts({inserts[0]});
    parlay::sequence<node_id> remaining_inserts(inserts.size() - 1);

    parlay::parallel_for(0, inserts.size(), [&](size_t i) {
      if (i == 0) {
        init_inserts[i] = inserts[i];
      } else {
        remaining_inserts[i - 1] = inserts[i];
      }
    });

    std::cout << "Performing functional update: update size = "
              << init_inserts.size() << std::endl;
    Graph new_G = insert_functional(S.graph, n, init_inserts);
    std::cout << "Performing inplace update on node: "
              << ((size_t)new_G.get_root()) << std::endl;
    std::cout << "Performing inplace update of size: " << remaining_inserts.size() << std::endl;
    insert_inplace(new_G, n, remaining_inserts);
    modify_edges(new_G, prev_inserts, n);
    std::cout << "Posting new version " << std::endl;
    VG.add_version_from_graph(std::move(new_G));
    VG.release_version(std::move(S));
  }

  Graph insert_functional(Graph &G, node_id n, parlay::sequence<node_id> to_insert){
    //generate 10 random edges for each point (size at most n) and then batch insert
    parlay::random_generator gen;
    std::uniform_int_distribution<int> dis(0, n-1);
    auto KVs = parlay::tabulate(to_insert.size(), [&] (size_t i){
      auto r = gen[i];
      parlay::sequence<node_id> new_nbh(10);
      for(int j=0; j<10; j++){
        new_nbh[j] = dis(r);
      }
      node_id index = to_insert[i];
      auto new_nghs = parlay::make_slice(new_nbh.begin(), new_nbh.end());
      size_t nghs_size = new_nbh.size();
      auto begin =
          (std::tuple<node_id, empty_weight>*)(new_nbh.begin());
      auto tree = edge_tree(begin, begin + nghs_size);
      auto tree_ptr = tree.root;
      tree.root = nullptr;
      return std::make_tuple(index, tree_ptr);
    });
    Graph new_G = G.insert_vertices_batch_functional(KVs.size(), KVs.begin());
    return new_G;
  }

  void access_edges(Graph &G, parlay::sequence<node_id> to_access){
    parlay::parallel_for(0, to_access.size(), [&] (size_t i){
      auto vtx = G.get_vertex(to_access[i]);
      auto f = [&] (node_id u, node_id v, empty_weight wgh) {
        auto hh = v; 
        return true;
      };
      vtx.out_neighbors().foreach_cond(f);
    });
  }

  void modify_edges(Graph &G, parlay::sequence<node_id> to_modify, node_id n){
    parlay::random_generator gen;
    std::uniform_int_distribution<int> dis(0, n-1);
    auto KVs = parlay::tabulate(to_modify.size(), [&] (size_t i){
      auto r = gen[i];
      parlay::sequence<node_id> new_nbh;
      auto vtx = G.get_vertex(to_modify[i]);
      auto f = [&] (node_id u, node_id v, empty_weight wgh) {
        new_nbh.push_back(v);
        return true;
      };
      vtx.out_neighbors().foreach_cond(f);
      
      new_nbh.push_back(dis(r));

      node_id index = to_modify[i];
      auto new_nghs = parlay::make_slice(new_nbh.begin(), new_nbh.end());
      size_t nghs_size = new_nbh.size();
      auto begin =
          (std::tuple<node_id, empty_weight>*)(new_nbh.begin());
      auto tree = edge_tree(begin, begin + nghs_size);
      auto tree_ptr = tree.root;
      tree.root = nullptr;
      return std::make_tuple(index, tree_ptr);
    });
    G.insert_vertices_batch(KVs.size(), KVs.begin());
  }

  void insert_inplace(Graph &G, node_id n, parlay::sequence<node_id> to_insert){
    parlay::random_generator gen;
    std::uniform_int_distribution<int> dis(0, n-1);
    auto KVs = parlay::tabulate(to_insert.size(), [&] (size_t i){
      auto r = gen[i];
      parlay::sequence<node_id> new_nbh(10);
      for(int j=0; j<10; j++){
        new_nbh[j] = dis(r);
      }
      node_id index = to_insert[i];
      auto new_nghs = parlay::make_slice(new_nbh.begin(), new_nbh.end());
      size_t nghs_size = new_nbh.size();
      auto begin =
          (std::tuple<node_id, empty_weight>*)(new_nbh.begin());
      auto tree = edge_tree(begin, begin + nghs_size);
      auto tree_ptr = tree.root;
      tree.root = nullptr;
      return std::make_tuple(index, tree_ptr);
    });
    G.insert_vertices_batch(KVs.size(), KVs.begin());
  }

  void build(node_id p){
    Graph Initial_Graph;
    Initial_Graph.insert_vertex_inplace(p, nullptr);
    VG = aspen::versioned_graph<Graph>(std::move(Initial_Graph));
  }

  
  void query(parlay::sequence<node_id> to_query){
    auto S = VG.acquire_version();
    access_edges(S.graph, to_query);
    VG.release_version(std::move(S));
  }
};

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int maxDeg, int beamSize,
         double alpha, Distance* D) {
  parlay::internal::timer t("ANN", report_stats);
  {
    node_id max_size = 1000000;
    
    test_graph TG;
    TG.build(0);
    auto inserts = parlay::tabulate(100000, [&] (node_id j) {return j+1;});
    auto prev_edges = {(node_id) 0};
    TG.insert(max_size, inserts, prev_edges);

    auto updater = [&] () {
      for(int i=1; i<10; i++){
        auto prev_inserts = parlay::tabulate(100000, [&] (node_id j) {return j;});
        auto inserts = parlay::tabulate(100000, [&] (node_id j) {return 100000* (node_id) (i)+1+j;});
        TG.insert(max_size, inserts, prev_inserts);
      }
    };

    auto queries = [&] () {
      for(int i=0; i<30; i++){
        auto inserts = parlay::tabulate(100000, [&] (node_id j) {return j;});
        TG.query(inserts);
      }
    };

    parlay::par_do(updater, queries);
  };
}

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
  int beamSize, int Q, double alpha,
  parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<ivec_point> groundTruth, 
  char* res_file, Distance* D) {
  parlay::internal::timer t("ANN", report_stats);
  {
    

  };
}
