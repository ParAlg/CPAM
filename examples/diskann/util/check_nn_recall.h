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
// #include "indexTools.h"
#include <set>
#include "../types.h"
// #include "beamSearch.h"
#include "csvfile.h"
#include "parse_results.h"
// #include "index.h"

template<typename T>
nn_result checkRecall(
        parlay::sequence<Tvec_point<T>*> &q,
        parlay::sequence<ivec_point> groundTruth,
        knn_index<T> &I,
        int k,
        int beamQ,
        float cut) {
  parlay::internal::timer t;
  int r = 10;

  float query_time;
  parlay::sequence<parlay::sequence<node_id>> query_results(q.size());
  parlay::parallel_for(0, q.size(), [&](size_t i){query_results[i]=parlay::sequence<node_id>(k);});
  I.query(q, k, beamQ, cut, query_results);
  // std::cout << "first set" << std::endl;

  t.next_time();
  I.query(q, k, beamQ, cut, query_results);

  query_time = t.next_time();
  // std::cout << "second set" << std::endl;
  float recall = 0.0;
  // bool dists_present = (groundTruth[0].distances.size() != 0);
  // if (groundTruth.size() > 0 && !dists_present) {
    size_t n = q.size();
    int numCorrect = 0;
    for(size_t i=0; i<n; i++){
      std::set<int> reported_nbhs;
      for(int l=0; l<r; l++) reported_nbhs.insert(query_results[i][l]);
      for(int l=0; l<r; l++){
	      if (reported_nbhs.find((groundTruth[i].coordinates)[l]) != reported_nbhs.end()){
          numCorrect += 1;
      }
      }
    }
    recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
  // }else if(groundTruth.size() > 0 && dists_present){
  //   size_t n = q.size();
  //   int numCorrect = 0;
  //   for(int i=0; i<n; i++){
  //     parlay::sequence<int> results_with_ties;
  //     for(int l=0; l<r; l++) results_with_ties.push_back(groundTruth[i].coordinates[l]);
  //     float last_dist = groundTruth[i].distances[r-1];
  //     for(int l=r; l<groundTruth[i].coordinates.size(); l++){
  //       if(groundTruth[i].distances[l] == last_dist){ 
  //         results_with_ties.push_back(groundTruth[i].coordinates[l]);
  //       }
  //     }
  //     std::set<int> reported_nbhs;
  //     // for(int l=0; l<r; l++) reported_nbhs.insert((q[i]->ngh)[l]);
  //     std::cout << results_with_ties.size() << std::endl;
  //     for(int l=0; l<results_with_ties.size(); l++){
	//       if (reported_nbhs.find(results_with_ties[l]) != reported_nbhs.end()){
  //         numCorrect += 1;
  //     }
  //     }
  //   }
  //   recall = static_cast<float>(numCorrect)/static_cast<float>(r*n);
  // }
  float QPS = q.size()/query_time;
  // std::cout << "QPS: " << QPS << std::endl;
  // std::cout << "Recall: " << recall << std::endl;
  // auto stats = query_stats(q);
  nn_result N(recall, QPS, k, beamQ, cut, q.size());
  return N;
}

// void write_to_csv(std::string csv_filename, parlay::sequence<float> buckets, 
//   parlay::sequence<nn_result> results, Graph G){
//   csvfile csv(csv_filename);
//   csv << "GRAPH" << "Parameters" << "Size" << "Build time" << "Avg degree" << "Max degree" << endrow;
//   csv << G.name << G.params << G.size << G.time << G.avg_deg << G.max_deg << endrow;
//   csv << endrow;
//   csv << "Num queries" << "Target recall" << "Actual recall" << "QPS" << "Average Cmps" << 
//     "Tail Cmps" << "Average Visited" << "Tail Visited" << "k" << "Q" << "cut" << endrow;
//   for(int i=0; i<results.size(); i++){
//     nn_result N = results[i];
//     csv << N.num_queries << buckets[i] << N.recall << N.QPS << N.avg_cmps << N.tail_cmps <<  
//       N.avg_visited << N.tail_visited << N.k << N.beamQ << N.cut << endrow;
//   }
//   csv << endrow;
//   csv << endrow;
// }

// parlay::sequence<int> calculate_limits(size_t avg_visited){
//   parlay::sequence<int> L(9);
//   for(float i=1; i<10; i++){
//     L[i-1] = (int) (i *((float) avg_visited) * .1);
//   }
//   auto limits = parlay::remove_duplicates(L);
//   return limits;
// }    

template<typename T>
void search_and_parse(parlay::sequence<Tvec_point<T>*> &q, 
    parlay::sequence<ivec_point> groundTruth, knn_index<T> &I){
    // unsigned d = v[0]->coordinates.size();

    parlay::sequence<nn_result> results;
    std::vector<int> beams = {15, 20, 30, 50, 75, 100, 125, 250, 500};
    std::vector<int> allk = {10, 15, 20, 30, 50, 100};
    std::vector<float> cuts = {1.1, 1.125, 1.15, 1.175, 1.2, 1.25};
    for (float cut : cuts)
      for (float Q : beams) 
        results.push_back(checkRecall(q, groundTruth, I, 10, Q, cut));

    // std::cout << "here1" << std::endl;
    for (float cut : cuts)
      for (int kk : allk)
        results.push_back(checkRecall(q, groundTruth, I, kk, 500, cut));
    // std::cout << "here2" << std::endl;
    // // check "limited accuracy"
    // parlay::sequence<int> limits = calculate_limits(results[0].avg_visited);
    // for(int l : limits){
    //   results.push_back(checkRecall(v, q, groundTruth, 10, 15, 1.14, d, random, l, start_point, D));
    // }

    // check "best accuracy"
    results.push_back(checkRecall(q, groundTruth, I, 100, 1000, 10.0));
    // std::cout << "here3" << std::endl;
    parlay::sequence<float> buckets = {.1, .15, .2, .25, .3, .35, .4, .45, .5, .55, .6, .65, .7, .73, .75, .77, .8, .83, .85, .87, .9, .93, .95, .97, .99, .995, .999};
    auto [res, ret_buckets] = parse_result(results, buckets);
    // if(res_file != NULL) write_to_csv(std::string(res_file), ret_buckets, res, G);
}

