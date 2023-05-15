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

#include "../IO.h"
// #include "../index.h"
#include "../types.h"
<<<<<<< HEAD
// #include "../writeSeq.h"
=======
>>>>>>> 6b201d31f1299636a8119719c338df5e3ca66051
// #include "util/check_nn_recall.h"

bool report_stats;

<<<<<<< HEAD
// using namespace benchIO;

=======
>>>>>>> 6b201d31f1299636a8119719c338df5e3ca66051
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

// *************************************************************
//  TIMING
// *************************************************************

template<typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>> &pts,
  int rounds, int R, int beamSize, double alpha, Distance* D)
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});

  time_loop(rounds, 0,
  [&] () {},
  [&] () {
    ANN<T>(v, R, beamSize, alpha, D);
  },
  [&] () {});
}

template<typename T>
void timeNeighbors(parlay::sequence<Tvec_point<T>> &pts,
		   parlay::sequence<Tvec_point<T>> &qpoints,
		   int k, int rounds, int R, int beamSize,
		   int beamSizeQ, double alpha,
		   parlay::sequence<ivec_point>& groundTruth, char* res_file, Distance* D)
{
  size_t n = pts.size();
  auto v = parlay::tabulate(n, [&] (size_t i) -> Tvec_point<T>* {
      return &pts[i];});

  size_t q = qpoints.size();
  auto qpts =  parlay::tabulate(q, [&] (size_t i) -> Tvec_point<T>* {
      return &qpoints[i];});

    time_loop(rounds, 0,
      [&] () {},
      [&] () {
        ANN<T>(v, k, R, beamSize, beamSizeQ, alpha, qpts, groundTruth, res_file, D);
      },
      [&] () {});

}

// Infile is a file in .fvecs format
int main(int argc, char* argv[]) {
    commandLine P(argc,argv,
    "[-a <alpha>] [-R <deg>]"
        "[-L <bm>] [-k <k> ] [-Q <bmq>] [-q <qF>]"
        "[-res <rF>] [-r <rnds>] [-t <tp>] [-D <df>] <inFile>");

  char* iFile = P.getArgument(0);
  char* qFile = P.getOptionValue("-q");
  char* cFile = P.getOptionValue("-c");
  char* rFile = P.getOptionValue("-res");
  char* vectype = P.getOptionValue("-t");
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
  char* dfc = P.getOptionValue("-D");

  distance_calls.reset();

  std::string df = std::string(dfc);
  Distance* D;
  if(df == "Euclidian") D = new Euclidian_Distance();
  else if(df == "mips") D = new Mips_Distance();
  else{
    std::cout << "Error: invalid distance type" << std::endl;
    abort();
  }

  std::string tp = std::string(vectype);

  if((tp != "uint8") && (tp != "int8") && (tp != "float")){
    std::cout << "Error: vector type not specified correctly, specify int8, uint8, or float" << std::endl;
    abort();
  }

  parlay::sequence<ivec_point> groundTruth;
  if(cFile != NULL) groundTruth = parse_ibin(cFile);
  
  if(tp == "float"){
    auto points = parse_fbin(iFile);
    if(qFile != NULL){
      auto qpoints = parse_fbin(qFile);
      timeNeighbors<float>(points, qpoints, k, rounds, R, L, Q,
        alpha, groundTruth, rFile, D);
    }
    else timeNeighbors<float>(points, rounds, R, L, alpha, D);
  } else if(tp == "uint8"){
    auto points = parse_uint8bin(iFile);
    if(qFile != NULL){
      auto qpoints = parse_uint8bin(qFile);
      timeNeighbors<uint8_t>(points, qpoints, k, rounds, R, L, Q,
        alpha, groundTruth, rFile, D);
    }
    else timeNeighbors<uint8_t>(points, rounds, R, L, alpha, D);
  } else if(tp == "int8"){
    auto points = parse_int8bin(iFile);
    if(qFile != NULL){
      auto qpoints = parse_int8bin(qFile);
      timeNeighbors<int8_t>(points, qpoints, k, rounds, R, L, Q,
        alpha, groundTruth, rFile, D);
    }
    else timeNeighbors<int8_t>(points, rounds, R, L, alpha, D);
  }
  
}