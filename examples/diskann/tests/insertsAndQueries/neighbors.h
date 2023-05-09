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
    unsigned d = (v[0]->coordinates).size();
    std::cout << "Size of dataset: " << v.size() << std::endl;
    using findex = knn_index<T>;
    findex I(v, maxDeg, beamSize, alpha, d, D);
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

    int k=10;
    int Q=100;
    auto queries = [&] () {
        // timer query_t;
        for(int i=0; i< (int) num_query_batches; i++){
            std::cout << "Querying elements " << query_start+(i*query_batch_size) 
            << " through " << query_start+((i+1)*query_batch_size)  << std::endl;
            auto queries = parlay::tabulate(query_batch_size, [&] (size_t j){
                return v[query_start + i*query_batch_size+j];});
            I.query(queries, k, Q);
            std::cout << "Finished query batch" << std::endl;
        }
    };

    auto empty = [&] () {};

    parlay::par_do(updater, queries);

    // parlay::par_do(updater, empty);

    // parlay::par_do(queries, empty);
}

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*> &v, int k, int maxDeg,
  int beamSize, int Q, double alpha,
  parlay::sequence<Tvec_point<T>*> &q, parlay::sequence<ivec_point> groundTruth, 
  char* res_file, Distance* D) {}

