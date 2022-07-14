#include <cstring>

#include <vector>
#include <algorithm>
#include <chrono>
#include <thread>
#include <cmath>

#include <iostream>
#include <fstream>

#include "aspen/aspen.h"
#include <cpam/parse_command_line.h>
#include <parlay/random.h>
#include "pbbsrandom.h"
#include "rmat_util.h"

using namespace std;

namespace aspen {

using namespace parlay;

template <class Graph>
double BatchUpdates_runner(Graph& G, cpam::commandLine P) {
  size_t n = G.num_vertices();

  auto update_sizes = parlay::sequence<size_t>(10);
  update_sizes[0] = 10;
  update_sizes[1] = 100;
  update_sizes[2] = 1000;
  update_sizes[3] = 10000;
  update_sizes[4] = 100000;
  update_sizes[5] = 1000000;
  update_sizes[6] = 10000000;
  update_sizes[7] = 100000000;
  update_sizes[8] = 1000000000;
  update_sizes[9] = 2000000000;

//  auto update_sizes = parlay::sequence<size_t>(3);
//  update_sizes[0] = 100000000;
//  update_sizes[1] = 1000000000;
//  update_sizes[2] = 2000000000;


  auto update_times = std::vector<double>();
  size_t n_trials = 9;
  pbbs::random r;

#ifdef USE_PAM_UPPER
  using pair_vertex = std::pair<uintE, uintE>;
#else
  using pair_vertex = std::tuple<uintE, uintE>;
#endif

  size_t start = 0;
  for (size_t us=start; us<update_sizes.size(); us++) {
    std::vector<double> inserts;
    std::vector<double> deletes;
    cout << "Running bs: " << update_sizes[us] << endl;

    for (size_t ts=0; ts<n_trials; ts++) {
      size_t updates_to_run = update_sizes[us];
      auto updates = parlay::sequence<pair_vertex>(updates_to_run);

      double a = 0.5;
      double b = 0.1;
      double c = 0.1;
      size_t nn = 1 << (parlay::log2_up(n) - 1);
      auto rmat = rMat<uintE>(nn, r.ith_rand(0), a, b, c);

      parallel_for(0, updates.size(), [&] (size_t i) {
        updates[i] = rmat(i);
      });

      {
        timer st; st.start();
        //G.insert_edges_batch(update_sizes[us], updates.begin(), false, true, nn, false);
        G.insert_edges_batch_inplace(update_sizes[us], updates.begin(), false, true, nn, false);
        //G = G.insert_edges_batch(update_sizes[us], updates.begin(), false, true, nn, false);
        double batch_time = st.stop();

        inserts.push_back(batch_time);
      }

      {
        timer st; st.start();
        //G.delete_edges_batch(update_sizes[us], updates.begin(), false, true, nn, false);
        G.delete_edges_batch_inplace(update_sizes[us], updates.begin(), false, true, nn, false);
        //G = G.delete_edges_batch(update_sizes[us], updates.begin(), false, true, nn, false);
        double batch_time = st.stop();

        deletes.push_back(batch_time);
      }
    }

    std::sort(inserts.begin(), inserts.end());
    std::sort(deletes.begin(), deletes.end());

    cout << "RESULT: Insert," << update_sizes[us] << "," << inserts[n_trials/2] << endl;
    cout << "RESULT: Delete," << update_sizes[us] << "," << deletes[n_trials/2] << endl;
  }
  return 1.0;
}


}  // namespace aspen

generate_symmetric_aspen_main(aspen::BatchUpdates_runner, false);
