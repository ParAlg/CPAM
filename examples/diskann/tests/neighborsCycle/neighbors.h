#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

#include <parlay/io.h>
#include <parlay/primitives.h>
#include <parlay/random.h>

#include <pam/get_time.h>
#include <pam/parse_command_line.h>

#include "../../index.h"
#include "../../util/check_nn_recall.h"


template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int maxDeg, int beamSize,
         double alpha, Distance* D) {
  parlay::internal::timer t("ANN", report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d, D);
    timer build_t;
    build_t.start();
    size_t n = v.size();
    I.build_index(parlay::tabulate(
        n/2, [&](size_t i) { return static_cast<node_id>(i); }));
    build_t.next("Build time for half the points");
    I.print_graph_status();
    auto inserts = parlay::tabulate(
        n/2, [&](size_t i) { return static_cast<node_id>(i+n/2); });
    I.insert(inserts);
    I.print_graph_status();
    build_t.next("Finished insertions");
  };
}



template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
  int beamSize, int Q, double alpha,
  parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<ivec_point> groundTruth, 
  char* res_file, Distance* D) {
  parlay::internal::timer t("ANN", report_stats);
  {
    timer build_t;
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d, D);
    build_t.start();
    // size_t n = v.size();
    I.build_index(parlay::tabulate(
        v.size(), [&](size_t i) { return static_cast<node_id>(i); }));
    build_t.next("Build time");
    int parts = 20;
    size_t m = v.size()/parts;
    for(int i=0; i<parts; i++){
      parlay::sequence<node_id> indices = parlay::tabulate(m, [&] (size_t j){return static_cast<node_id>(i*m+j);});
      std::cout << "Deleting indices " << indices[0] << " through " << indices[m-1] << std::endl; 
      I.lazy_delete(indices);
      I.start_delete_epoch();
      I.consolidate_deletes(parlay::tabulate(
        v.size(), [&](size_t i) { return static_cast<node_id>(i); }));
      I.end_delete_epoch();
      std::cout << "Re-inserting " << indices[0] << " through " << indices[m-1] << std::endl; 
      I.insert(indices);
    }
    build_t.next("Finished rebuilding");
    // std::cout << "Now performing queries" << std::endl;
    search_and_parse(q, groundTruth, I);
    build_t.next("query time");

   
  };
}

