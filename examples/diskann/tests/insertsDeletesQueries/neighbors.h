#include <pam/get_time.h>
#include <pam/parse_command_line.h>
#include <parlay/io.h>
#include <parlay/primitives.h>
#include <parlay/random.h>

#include <algorithm>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <set>
#include <vector>

#include "../../index.h"
#include "../../util/check_nn_recall.h"

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*>& v, int maxDeg, int beamSize,
         double alpha, Distance* D) {
  unsigned d = (v[0]->coordinates).size();

  using findex = knn_index<T>;
  findex I(v, maxDeg, beamSize, alpha, d, D);

  size_t n = v.size();
  size_t insert_batch_size = 10000;
  size_t query_batch_size = 20000;

  float build_frac = .6;
  float delete_frac = .1;
  float insert_frac = .2;
  float query_frac = .2;
  

  size_t num_inserts = n * insert_frac;
  size_t num_queries = n * query_frac;
  size_t num_deletes = n * delete_frac;
  size_t num_build = n * build_frac;

  size_t num_insert_batches = num_inserts / insert_batch_size;
  size_t num_query_batches = num_queries / query_batch_size;

  size_t delete_batch_size = num_build / num_insert_batches;

  size_t query_start = n * (build_frac + insert_frac);
  size_t insert_start = n * (build_frac);
  node_id delete_start = n* (build_frac - delete_frac);
  // std::cout << query_start << std::endl;

  auto build_indices = parlay::tabulate(num_build, [&] (node_id i){return i;});
  I.build_index(build_indices);

  int k=10;
  int Q=100;

  auto updater = [&]() {
    // timer update_t;
    auto to_delete = parlay::tabulate(num_deletes, [&] (node_id i){return delete_start+i;});
    I.lazy_delete(to_delete);
    I.start_delete_epoch();
    
    for (size_t i = 0; i < num_insert_batches; i++) {
      parlay::sequence<node_id> indices = parlay::tabulate(insert_batch_size, [&](node_id j) {
          return static_cast<node_id>(insert_start + i * insert_batch_size + j);
      });
      std::cout << "Inserting indices " << indices[0] << " through "
                << indices[indices.size() - 1] << std::endl;
      I.insert(indices);
      std::cout << "Finished inserting" << std::endl;
      parlay::sequence<node_id> to_consolidate = parlay::tabulate(delete_batch_size, [&](node_id j) {
          return static_cast<node_id>(i * delete_batch_size + j);
      });
      std::cout << "Consolidating for indices " << to_consolidate[0] << " through " << to_consolidate[to_consolidate.size()-1] << std::endl;
      I.consolidate_deletes(to_consolidate);
      std::cout << "Finished consolidating" << std::endl;
    }
    I.end_delete_epoch();
  };

  auto queries = [&]() {
    // timer query_t;
    for (int i = 0; i < (int)num_query_batches; i++) {
      std::cout << "Querying elements " << query_start + (i * query_batch_size)
                << " through " << query_start + ((i + 1) * query_batch_size)
                << std::endl;
      auto queries = parlay::tabulate(query_batch_size, [&](size_t j) {
        return v[query_start + i * query_batch_size + j];
      });
      I.query(queries, k, Q);
      std::cout << "Finished query batch" << std::endl;
      std::this_thread::sleep_for(std::chrono::milliseconds(1000));
    }
  };

  auto empty = [&]() {};


  parlay::par_do(updater, queries);

  // parlay::par_do(updater, empty);

  // parlay::par_do(queries, empty);
}

template <typename T>
void ANN(parlay::sequence<Tvec_point<T>*>& v, int k, int maxDeg, int beamSize,
         int Q, double alpha, parlay::sequence<Tvec_point<T>*>& q,
         parlay::sequence<ivec_point> groundTruth, char* res_file,
         Distance* D) {
  unsigned d = (v[0]->coordinates).size();

  using findex = knn_index<T>;
  findex I(v, maxDeg, beamSize, alpha, d, D);

  size_t n = v.size();
  size_t insert_batch_size = 10000;
  size_t query_batch_size = 20000;

  float build_frac = .6;
  float delete_frac = .1;
  float insert_frac = .2;
  float query_frac = .2;
  

  size_t num_inserts = n * insert_frac;
  size_t num_queries = n * query_frac;
  size_t num_deletes = n * delete_frac;
  size_t num_build = n * build_frac;

  size_t num_insert_batches = num_inserts / insert_batch_size;
  size_t num_query_batches = num_queries / query_batch_size;

  size_t delete_batch_size = num_build / num_insert_batches;

  size_t query_start = n * (build_frac + insert_frac);
  size_t insert_start = n * (build_frac);
  node_id delete_start = n* (build_frac - delete_frac);
  // std::cout << query_start << std::endl;

  auto build_indices = parlay::tabulate(num_build, [&] (node_id i){return i;});
  I.build_index(build_indices);

//   auto updater = [&]() {
//     // timer update_t;
//     auto to_delete = parlay::tabulate(num_deletes, [&] (node_id i){return delete_start+i;});
//     I.lazy_delete(to_delete);
//     I.start_delete_epoch();
    
//     for (size_t i = 0; i < num_insert_batches; i++) {
//       parlay::sequence<node_id> indices = parlay::tabulate(insert_batch_size, [&](node_id j) {
//           return static_cast<node_id>(insert_start + i * insert_batch_size + j);
//       });
//       std::cout << "Inserting indices " << indices[0] << " through "
//                 << indices[indices.size() - 1] << std::endl;
//       I.insert(indices);
//       std::cout << "Finished inserting" << std::endl;
//       parlay::sequence<node_id> to_consolidate = parlay::tabulate(delete_batch_size, [&](node_id j) {
//           return static_cast<node_id>(delete_start + i * delete_batch_size + j);
//       });
//       std::cout << "Consolidating for indices " << to_consolidate[0] << " through " << to_consolidate[to_consolidate.size()-1] << std::endl;
//       I.consolidate_deletes(to_consolidate);
//       std::cout << "Finished consolidating" << std::endl;
//     }
//     I.end_delete_epoch();
//   };

//   auto queries = [&]() {
//     // timer query_t;
//     for (int i = 0; i < (int)num_query_batches; i++) {
//       std::cout << "Querying elements " << query_start + (i * query_batch_size)
//                 << " through " << query_start + ((i + 1) * query_batch_size)
//                 << std::endl;
//       auto queries = parlay::tabulate(query_batch_size, [&](size_t j) {
//         return v[query_start + i * query_batch_size + j];
//       });
//       I.query(queries, k, Q);
//       std::cout << "Finished query batch" << std::endl;
//       std::this_thread::sleep_for(std::chrono::milliseconds(1000));
//     }
//   };

//   auto empty = [&]() {};

// //  updater();
// //
// //  std::cout << "Finished updates" << std::endl;
// //
// //  queries();

//   parlay::par_do(updater, queries);

  // parlay::par_do(updater, empty);

  // parlay::par_do(queries, empty);

}
