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
    int parts = 10; 
    size_t n = v.size();
    size_t m = (size_t) (n/parts);
    size_t r = m/2; 
    parlay::sequence<node_id> inserts = parlay::tabulate(m, [&] (size_t i){return static_cast<node_id>(n-i-1);});
    std::cout << "Building with points " << inserts[inserts.size()-1] << " through " << inserts[0] << std::endl; 
    I.build_index(inserts);
    for(int i=parts*2-1; i>1; i--){
      parlay::sequence<node_id> delete_list = parlay::tabulate(r, [&] (size_t j){return static_cast<node_id>(((i)*r)+r-j-1);});
      std::cout << "Deleting points " << delete_list[delete_list.size()-1] << " through " << delete_list[0] << std::endl; 
      I.lazy_delete(delete_list);
      I.start_delete_epoch();
      parlay::sequence<node_id> current_list = parlay::tabulate(r, [&] (size_t j){
        return static_cast<node_id>((i-1)*(r)+j);});
      std::cout << "Consolidating based on ids " << current_list[0] << " to " << current_list[current_list.size()-1] << std::endl;
      I.consolidate_deletes(current_list);
      parlay::sequence<node_id> medoid_list = {I.get_medoid()};
      I.consolidate_deletes(medoid_list);
      I.end_delete_epoch();
      parlay::sequence<node_id> insert_list = parlay::tabulate(r, [&] (size_t j){return static_cast<node_id>(((i)*r)-r-j-1);});
      std::cout << "Inserting points " << insert_list[insert_list.size()-1] << " through " << insert_list[0] << std::endl; 
      I.insert(insert_list);
    }
  };
}

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
    int parts = 10; 
    size_t n = v.size();
    size_t m = (size_t) (n/parts);
    size_t r = m/2; 
    std::cout << "Here" << std::endl;
    parlay::sequence<node_id> inserts = parlay::tabulate(m, [&] (size_t i){return static_cast<node_id>(n-i-1);});
    std::cout << "Building with points " << inserts[inserts.size()-1] << " through " << inserts[0] << std::endl; 
    I.build_index(inserts);
    for(int i=parts*2-1; i>1; i--){
      parlay::sequence<node_id> delete_list = parlay::tabulate(r, [&] (size_t j){return static_cast<node_id>(((i)*r)+r-j-1);});
      std::cout << "Deleting points " << delete_list[delete_list.size()-1] << " through " << delete_list[0] << std::endl; 
      I.lazy_delete(delete_list);
      parlay::sequence<node_id> current_list = parlay::tabulate(r, [&] (size_t j){
        return static_cast<node_id>(i*r+j-1);});
      I.start_delete_epoch();
      I.consolidate_deletes(current_list);
      I.end_delete_epoch();
      parlay::sequence<node_id> insert_list = parlay::tabulate(r, [&] (size_t j){return static_cast<node_id>(((i)*r)-r-j-1);});
      std::cout << "Inserting points " << insert_list[insert_list.size()-1] << " through " << insert_list[0] << std::endl; 
      I.insert(insert_list);
    }

    parlay::sequence<parlay::sequence<unsigned>> query_results(q.size());
    std::cout << "Built index, now performing queries" << std::endl;
    search_and_parse(q, groundTruth, I);

  };
}
