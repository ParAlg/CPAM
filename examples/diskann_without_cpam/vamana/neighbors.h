// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parlay/random.h"
//#include "common/geometry.h"
#include "../utils/NSGDist.h"
#include "../utils/types.h"
#include "index.h"
#include "../utils/beamSearch.h"
#include "../utils/indexTools.h"
#include "../utils/stats.h"

extern bool report_stats;

template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg, int beamSize, int beamSizeQ, double alpha, double dummy,

  parlay::sequence<Tvec_point<T>*> &q) {
  parlay::internal::timer t("ANN",report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(maxDeg, beamSize, alpha, d);


    int parts = 10;
    size_t n = v.size();
    size_t m = (size_t) (n/parts);
    size_t r = m/2;
    parlay::sequence<int> inserts = parlay::tabulate(n, [&] (size_t i){return static_cast<int>(i);});
    std::cout << "Building with points " << inserts[inserts.size()-1] << " through " << inserts[0] << std::endl;
    I.build_index(v, inserts);

    I.searchNeighbors(q, v, beamSizeQ, k);
    t.next("Found nearest neighbors");
    if(report_stats){
      //average numbers of nodes searched using beam search
      graph_stats(v);
      query_stats(q);
      t.next("stats");
    }
  };
}


template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize, double alpha, double dummy) {
  parlay::internal::timer t("ANN",report_stats);
  distance_calls.reset();
  total_visited.reset();
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(maxDeg, beamSize, alpha, d);
    I.build_index(v, parlay::tabulate(v.size(), [&] (size_t i){return static_cast<int>(i);}));
    std::cout << "Total vertices visited: " << total_visited.get_value() << std::endl;
    std::cout << "Total distance calls: " << distance_calls.get_value() << std::endl;
    t.next("Built index");
    if(report_stats){
      graph_stats(v);
      t.next("stats");
    }
  };
}
