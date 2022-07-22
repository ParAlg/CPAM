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

bool report_stats;

template<typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize, double alpha, double dummy) {
  parlay::internal::timer t("ANN",report_stats);
  {
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;

    using findex = knn_index<T>;
    findex I(maxDeg, beamSize, alpha, d);
    I.build_index(v, parlay::tabulate(v.size(), [&] (size_t i){return static_cast<int>(i);}));
//    t.next("Built index");
//    if(report_stats){
//      graph_stats(v);
//      t.next("stats");
//    }
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
                   int beamSize, double delta, double alpha, char* outFile) {
  size_t n = pts.size();
  auto v =
      parlay::tabulate(n, [&](size_t i) -> Tvec_point<T>* { return &pts[i]; });

  time_loop(rounds, 0, [&]() {},
            [&]() { ANN<T>(v, maxDeg, beamSize, alpha, delta); }, [&]() {});
}

int main(int argc, char** argv) {
  commandLine P(argc, argv,
                "[-a <alpha>] [-d <delta>] [-R <deg>]"
                "[-L <bm>] [-k <k> ] [-Q <bmq>] [-q <qF>] [-o <oF>] [-r "
                "<rnds>] [-b <algoOpt>] <inFile>");

  char* iFile = P.getArgument(0);
  char* oFile = P.getOptionValue("-o");
  int R = P.getOptionIntValue("-R", 5);
  if (R < 1) P.badArgument();
  int L = P.getOptionIntValue("-L", 10);
  if (L < 1) P.badArgument();
  int Q = P.getOptionIntValue("-Q", L);
  if (Q < 1) P.badArgument();
  int rounds = P.getOptionIntValue("-r", 1);
  int k = P.getOptionIntValue("-k", 1);
  if (k > 1000 || k < 1) P.badArgument();
  double alpha = P.getOptionDoubleValue("-a", 1.2);
  double delta = P.getOptionDoubleValue("-d", .01);
  int HCNNG = P.getOptionIntValue("-b", 0);

  // TODO: add back qfile.

  bool fvecs = true;
  std::string filename = std::string(iFile);
  std::string::size_type n = filename.size();
  if (filename[n - 5] == 'b') fvecs = false;

  int maxDeg;
  if (HCNNG != 0)
    maxDeg = R * L;
  else
    maxDeg = R;

  if (fvecs) {  // vectors are floating point coordinates
    parlay::sequence<Tvec_point<float>> points = parse_fvecs(iFile, maxDeg);
    std::cout << "Parsed fvecs, len = " << points.size() << std::endl;
    timeNeighbors<float>(points, rounds, R, L, delta, alpha, oFile);
  } else {  // vectors are uint8 coordinates
    parlay::sequence<Tvec_point<uint8_t>> points = parse_bvecs(iFile, maxDeg);
    timeNeighbors<uint8_t>(points, rounds, R, L, delta, alpha, oFile);
  }
}
