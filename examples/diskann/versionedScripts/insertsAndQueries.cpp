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
#include "../index.h"
#include "../types.h"
#include "../writeSeq.h"

bool report_stats;

using namespace benchIO;

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> v, int maxDeg, int beamSize,
         double alpha, double dummy, int k, int Q) {

    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d);
    I.build_index({0});
    size_t n = v.size();
    size_t update_batch_size = 50000;
    size_t query_batch_size = 10000;

    float update_frac = .8;
    float query_frac = 1-update_frac;

    size_t num_updates = n*update_frac;
    size_t num_queries = n*query_frac;

    size_t num_update_batches = num_updates/update_batch_size;
    size_t num_query_batches = num_queries/query_batch_size;

    size_t query_start = n*update_frac;
    // std::cout << query_start << std::endl;

    auto updater = [&] () {
        // timer update_t;
        for(node_id i=0; i< (node_id) num_update_batches; i++){
            parlay::sequence<node_id> indices;
            if(i==0){ 
              indices = parlay::tabulate(update_batch_size-1, [&] (node_id j){
                return static_cast<node_id>(i*update_batch_size+j+1);});
            }else indices = parlay::tabulate(update_batch_size, [&] (node_id j){
                return static_cast<node_id>(i*update_batch_size+j);});
            std::cout << "Inserting indices " << indices[0] << 
                " through " << indices[indices.size()-1] << std::endl;
            I.insert(indices);
            std::cout << "Finished inserting" << std::endl;
        }
    };

    auto queries = [&] () {
        // timer query_t;
        // for(int i=0; i< (int) num_query_batches; i++){
        //     std::cout << "Querying elements" << query_start+(i*query_batch_size) 
        //     << " through " << query_start+((i+1)*query_batch_size)  << std::endl;
        //     auto queries = parlay::tabulate(query_batch_size, [&] (size_t j){
        //         return v[query_start + i*query_batch_size+j];});
        //     I.query(queries, k, Q);
        //     std::cout << "Finished query batch" << std::endl;
        // }
    };

    parlay::par_do(updater, queries);
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
                   int beamSizeQ) {
  size_t n = pts.size();
  auto v =
      parlay::tabulate(n, [&](size_t i) -> Tvec_point<T>* { return &pts[i]; });


time_loop(rounds, 0, [&]() {},
            [&]() {
            ANN<T>(v, maxDeg, beamSize, alpha, delta, k, beamSizeQ);
            },
            [&]() {});

}

int main(int argc, char** argv) {
  commandLine P(argc, argv,
                "[-a <alpha>] [-d <delta>] [-R <deg>]"
                "[-L <bm>] [-k <k> ] [-Q <bmq>] [-r "
                "<rnds>] <inFile>");

  char* iFile = P.getArgument(0);
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

  distance_calls.reset();

  bool fvecs = true;
  std::string filename = std::string(iFile);
  std::string::size_type n = filename.size();
  if (filename[n - 5] == 'b') fvecs = false;

  stats = P.getOptionValue("-stats");

  if (fvecs) {  // vectors are floating point coordinates
    parlay::sequence<Tvec_point<float>> points = parse_fvecs(iFile);
    std::cout << "Parsed fvecs, len = " << points.size() << std::endl;
    timeNeighbors<float>(points, rounds, R, L, alpha, delta, k, Q);
  } else {  // vectors are uint8 coordinates
    parlay::sequence<Tvec_point<uint8_t>> points = parse_bvecs(iFile);
    timeNeighbors<uint8_t>(points, rounds, R, L, alpha, delta, k, Q);
  }
}