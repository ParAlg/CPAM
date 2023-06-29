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
    I.build_index(parlay::tabulate(
        v.size(), [&](size_t i) { return static_cast<node_id>(i); }));
    build_t.next("Build time");
  };
}

//v, k, R, beamSize, beamSizeQ, alpha, qpts, groundTruth, res_file, D



template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
  int beamSize, int Q, double alpha,
  parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<ivec_point> groundTruth, 
  char* res_file, Distance* D) {
  parlay::internal::timer t("ANN", report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d, D);
    t.start();
    I.build_index(parlay::tabulate( v.size(), [&](size_t i) { return static_cast<node_id>(i); }));
    t.next("Build time");
    search_and_parse(q, groundTruth, I);
    t.next("query time");
  };
}

