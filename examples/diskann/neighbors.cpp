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

#include "IO.h"
#include "index.h"
#include "types.h"
#include "writeSeq.h"

bool report_stats;

using namespace benchIO;

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize,
         double alpha, double dummy) {
  parlay::internal::timer t("ANN", report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d);
    timer build_t;
    build_t.start();
    I.build_index(parlay::tabulate(
        v.size(), [&](size_t i) { return static_cast<node_id>(i); }));
    build_t.next("Build time");
    // int parts = 20;
    // size_t m = v.size()/parts;
    // for(int i=0; i<parts; i++){
    //   parlay::sequence<node_id> indices = parlay::tabulate(m, [&] (size_t j){return static_cast<node_id>(i*m+j);});
    //   I.lazy_delete(indices);
    //   I.consolidate_deletes();
    //   I.insert(indices);
    // }
  };
}

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize,
         double alpha, double dummy, int k, int Q,
         parlay::sequence<Tvec_point<T>*> q, char* outFile) {
  parlay::internal::timer t("ANN", report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d);
    I.build_index(parlay::tabulate( v.size(), [&](size_t i) { return static_cast<node_id>(i); }));

    int parts = 20;
    size_t m = v.size()/parts;
    for(int i=0; i<parts; i++){
      parlay::sequence<node_id> indices = parlay::tabulate(m, [&] (size_t j){return static_cast<node_id>(i*m+j);});
      I.lazy_delete(indices);
      I.consolidate_deletes();
      I.insert(indices);
    }

    parlay::sequence<parlay::sequence<unsigned>> query_results(q.size());
    std::cout << "Built index, now performing queries" << std::endl;
    parlay::parallel_for(0, q.size(), [&](size_t i) {
      query_results[i] = I.query(q[i]->coordinates.begin(), k, Q);
    });

    if (outFile != NULL) {
      size_t m = q.size() * (k + 1);
      parlay::sequence<unsigned> Pout(m);
      parlay::parallel_for(0, q.size(), [&](size_t i) {
        Pout[(k + 1) * i] = q[i]->id;
        for (int j = 0; j < (int)query_results[i].size(); j++)
          Pout[(k + 1) * i + j + 1] = query_results[i][j];
        for (int j = query_results[i].size(); j < k; k++)
          Pout[(k + 1) * i + j + 1] = 0;
      });
      writeIntSeqToFile(Pout, outFile);
    }
  };
}

template <class F, class G, class H>
void time_loop(int rounds, double delay, F initf, G runf, H endf) {
  parlay::internal::timer t;
  // run for delay seconds to "warm things up"
  // will skip if delay is zero
  while (t.total_time() < delay) {
    initf();
    runf();
    endf();
  }
  for (int i = 0; i < rounds; i++) {
    initf();
    t.start();
    runf();
    t.next("");
    endf();
  }
}

template <typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>>& pts, int rounds, int maxDeg,
                   int beamSize, double alpha, double delta, int k,
                   int beamSizeQ, parlay::sequence<Tvec_point<T>>& queries,
                   char* outFile) {
  size_t n = pts.size();
  auto v =
      parlay::tabulate(n, [&](size_t i) -> Tvec_point<T>* { return &pts[i]; });

  if (queries.size() != 0) {
    std::cout << "here" << std::endl;
    auto q = parlay::tabulate(queries.size(), [&](size_t i) -> Tvec_point<T>* {
      return &queries[i];
    });
    time_loop(rounds, 0, [&]() {},
              [&]() {
                ANN<T>(v, maxDeg, beamSize, alpha, delta, k, beamSizeQ, q,
                       outFile);
              },
              [&]() {});
  } else {
    time_loop(rounds, 0, [&]() {},
              [&]() { ANN<T>(v, maxDeg, beamSize, alpha, delta); }, [&]() {});
  }
}

int main(int argc, char** argv) {
  commandLine P(argc, argv,
                "[-a <alpha>] [-d <delta>] [-R <deg>]"
                "[-L <bm>] [-k <k> ] [-Q <bmq>] [-q <qF>] [-o <oF>] [-r "
                "<rnds>] <inFile>");

  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  char* qFile = P.getOptionValue("-q");
  int R = P.getOptionIntValue("-R", 15);
  if (R < 1) P.badArgument();
  int L = P.getOptionIntValue("-L", 25);
  if (L < 1) P.badArgument();
  int Q = P.getOptionIntValue("-Q", L);
  if (Q < 1) P.badArgument();
  int rounds = P.getOptionIntValue("-r", 1);
  int k = P.getOptionIntValue("-k", 1);
  if (k > 1000 || k < 1) P.badArgument();
  double alpha = P.getOptionDoubleValue("-a", 1.2);
  double delta = P.getOptionDoubleValue("-d", .01);

  // TODO: add back qfile.

  bool fvecs = true;
  std::string filename = std::string(iFile);
  std::string::size_type n = filename.size();
  if (filename[n - 5] == 'b') fvecs = false;

  stats = P.getOptionValue("-stats");

  if (fvecs) {  // vectors are floating point coordinates
    parlay::sequence<Tvec_point<float>> points = parse_fvecs(iFile);
    std::cout << "Parsed fvecs, len = " << points.size() << std::endl;
    parlay::sequence<Tvec_point<float>> queries;
    if (qFile != NULL) {
      queries = parse_fvecs(qFile);
      std::cout << "Parsed queries, len = " << queries.size() << std::endl;
    }
    timeNeighbors<float>(points, rounds, R, L, alpha, delta, k, Q, queries,
                         oFile);
  } else {  // vectors are uint8 coordinates
    parlay::sequence<Tvec_point<uint8_t>> points = parse_bvecs(iFile);
    parlay::sequence<Tvec_point<uint8_t>> queries;
    if (qFile != NULL) queries = parse_bvecs(qFile);
    timeNeighbors<uint8_t>(points, rounds, R, L, alpha, delta, k, Q, queries,
                           oFile);
  }
}
