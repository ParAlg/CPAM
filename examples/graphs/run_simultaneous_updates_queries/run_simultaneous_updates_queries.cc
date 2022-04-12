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
#include <random>
#include <parlay/random.h>
#include "pbbsrandom.h"
#include "rmat_util.h"

# include "../algorithms/BFS.h"
# include "../algorithms/BC.h"
# include "../algorithms/MIS.h"

using namespace std;

namespace aspen {

using namespace parlay;

#ifdef USE_PAM_UPPER
  using pair_vertex = std::pair<uintE, uintE>;
#else
  using pair_vertex = std::tuple<uintE, uintE>;
#endif

template <class Seq>
void write_datapoints_to_csv(std::string& fname, Seq const& S, size_t n_points) {

  if (S.size() == n_points) {
    ofstream outf;
    outf.open(fname);
    outf << "time, latency" << endl;
    for (size_t i=0; i<n_points; i++) {
      outf <<  fixed << setprecision(9) << S[i].first << "," << S[i].second << endl;
    }
    outf.close();
    return;
  }

  double min_t = S[0].first;
  double max_t = S[S.size()-1].first;

  cout << "min_t = " << min_t << " max_t = " << max_t << " nm data points = " << S.size() << endl;

  auto output_pts = parlay::sequence<pair<double, double>>(n_points);
  double window_size = (max_t - min_t) / n_points;
  std::cout << "window_size = " << window_size << std::endl;

  size_t k = 0;
  double window_start = min_t;
  double window_end = min_t + window_size;
  size_t pts_in_last_window = 0;
  double total_in_last_window = 0.0;

  std::vector<double> pts_in_window;

  for (size_t i=0; i<S.size(); i++) {
    double t_i = S[i].first;
    double lat_i = S[i].second;
    if (t_i < window_end) {
      pts_in_last_window++;
      total_in_last_window += lat_i;
      pts_in_window.push_back(lat_i);
    } else {
      if (pts_in_last_window) {
        double avg_latency = total_in_last_window / pts_in_last_window;
        output_pts[k++] = make_pair((window_start + window_end) / 2, avg_latency);
      }
      window_start = window_end;
      window_end += window_size;
      pts_in_window.clear();

      pts_in_window.push_back(lat_i);
    }
  }

  ofstream outf;
  outf.open(fname);
  outf << "time, latency" << endl;
  for (size_t i=0; i<k; i++) {
    outf <<  fixed << setprecision(9) << output_pts[i].first << "," << output_pts[i].second << endl;
  }
  outf.close();
}

template <class Graph>
double SimultaneousUpdatesQueries_runner(Graph& G, cpam::commandLine P) {
  size_t n = G.num_vertices();
  size_t gran = 10;

  G.set_num_vertices(131341);
  pbbs::random r;
  double a = 0.5;
  double b = 0.1;
  double c = 0.1;
  size_t nn = 1 << (parlay::log2_up(n) - 1);
  auto rmat = rMat<uintE>(nn, r.ith_rand(0), a, b, c);

  size_t updates_to_run = std::max(n, 1000000UL);
  std::cout << "n = " << n << " updates to run = " << updates_to_run << " nn = " << nn << std::endl;
  auto updates = parlay::sequence<pair_vertex>(2*updates_to_run);

  parlay::parallel_for(0, updates_to_run, [&] (size_t i) {
    auto [u,v] = rmat(i);
    updates[2*i] = pair_vertex(u,v);
    updates[2*i+1] = pair_vertex(v,u);
  });

  std::cout << "G ref = " << G.ref_cnt() << " root = " << G.get_root()  << " num_vertices = " << G.num_vertices() << std::endl;

  bool query_only = P.getOption("-query_only");
  bool update_only = P.getOption("-update_only");
  auto algo_name = P.getOptionValue("-alg", "BFS");
  auto root = G.get_root();
  auto VG = versioned_graph<Graph>(std::move(G));
  std::cout << "Initially, timestamp is: " << VG.latest_timestamp() << std::endl;

  std::cout << "After creating vg root ref_cnt = " << versioned_graph<Graph>::Tree::ref_cnt(root) << std::endl;

  std::cout << "num_workers = " << parlay::num_workers() << std::endl;
  std::cout << "worker_id = " << parlay::worker_id() << std::endl;

  string update_fname = P.getOptionValue("-update-file", "updates.dat");
  string query_fname = P.getOptionValue("-query-file", "queries.dat");

  auto times = parlay::sequence<pair<double, double>>(updates_to_run/gran);

  auto acquire_version = [&] () {
    std::cout << "about to acquire..." << std::endl;
    auto vv = VG.acquire_version();
    vv.graph.set_num_vertices(n);
    return vv;
  };

  auto v1 = acquire_version();
  std::cout << "v1, n = " << v1.graph.num_vertices() << " root = " << v1.graph.get_root() << " ref_cnt = " << v1.graph.ref_cnt() << std::endl;
  std::cout << "v1, n = " << v1.graph.num_vertices() << std::endl;

  bool updates_done = false;
  double update_time;

  auto updater = [&] () {
    timer t("Updates");
    double last = 0.0;
    for (size_t i=0; i<(updates_to_run/gran); i++) {
      if (i % 10000 == 0) {
        auto vv = acquire_version();
        std::cout << "Finished " << (i*gran) << " updates. Num edges in graph is: " << vv.graph.num_edges() << std::endl;
        VG.release_version(std::move(vv));
      }
      auto slice = updates.cut(gran*i, gran*i+gran);
      VG.insert_edges_batch(slice);
      double elapsed = t.get_total();
      times[i] = std::make_pair(elapsed, elapsed - last);
      last = elapsed;
    }
    update_time = t.stop();
    updates_done = true;
  };

  std::default_random_engine generator;
  std::poisson_distribution<int> distribution(300); // mean is 300ms
  size_t queries_run = 0;
  double total_query_time = 0;

  std::vector<std::pair<double, double>> query_times;

  size_t bfs_run = 0;
  auto queries = [&] () {
    timer t; t.start();
    size_t max_queries = 100;
    size_t queries_done = 0;
    while (!updates_done || (query_only && queries_run < max_queries)) {
      int poisson_sleep = distribution(generator);
      auto S = acquire_version();

      timer bt; bt.start();
      std::cout << "Starting BFS on graph, ref_cnt = " << S.graph.ref_cnt() << " num_edges = " << S.graph.num_edges() << std::endl;
      if (algo_name == "BFS") {
        BFS(S.graph, 10012, /*flatsnap=*/true);
      } else if (algo_name == "BC") {
        BC(S.graph, 10012, /*flatsnap=*/true);
      } else {
        MaximalIndependentSet_rootset::MaximalIndependentSet(S.graph, /*flatsnap=*/true);
      }
      double tt = bt.stop();
      double elapsed = t.get_total();
      std::cout << "BFS took: " << tt << " elapsed = " << elapsed << std::endl;
      bfs_run++;

      query_times.emplace_back(elapsed, tt);

      VG.release_version(std::move(S));
      queries_run++;
      total_query_time += tt;

      std::this_thread::sleep_for(std::chrono::milliseconds(poisson_sleep));
      queries_done++;
      if (query_only && queries_done == max_queries) break;
    }
  };

 if (update_only) {
    std::cout << "Running update only" << std::endl;
    updater();
  } else if (query_only) {
    std::cout << "Running query only" << std::endl;
    queries();
  } else {
    std::cout << "Running concurrent updates + queries" << std::endl;
    parlay::par_do(updater, queries);
  }
  std::cout << "Finished " << queries_run << " many queries." << std::endl;

  std::cout << "At end of updates" << std::endl;
  std::cout << "v1 ref = " << v1.graph.ref_cnt() << " root = " << v1.graph.get_root() << std::endl;

  if (!query_only) {
    write_datapoints_to_csv(update_fname, times, 100);
  }
  if (!update_only) {
    write_datapoints_to_csv(query_fname, query_times, query_times.size());
  }
  if (!query_only) {
    std::cout << "RESULT: Update throughput = " << (updates_to_run / update_time) << std::endl;
  }
  if (!update_only) {
    std::cout << "RESULT: Average query time = " << (total_query_time / queries_run) << std::endl;
  }

  return 1.0;
}


}  // namespace aspen

int main(int argc, char* argv[]) {
  cpam::commandLine P(argc, argv, " [-s] <inFile>");

  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool mmap = P.getOptionValue("-m");
  if (!symmetric) {
    std::cout
        << "# The application expects the input graph to be symmetric (-s "
           "flag)."
        << std::endl;
    std::cout << "# Please run on a symmetric input." << std::endl;
  }
  timer rt;
  rt.start();
  std::cout << "About to start parsing graph!" << std::endl;
  auto G = aspen::parse_unweighted_symmetric_graph(iFile, mmap);
  rt.next("Graph read time");
  auto AG = aspen::symmetric_graph_from_static_graph(G);
  aspen::SimultaneousUpdatesQueries_runner(AG, P);
}
