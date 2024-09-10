#pragma once

#include "byte-pd-amortized.h"
#include "graph_io.h"
#include "immutable_graph.h"
#include "traversable_graph.h"
#include "utils.h"

#include <cpam/get_time.h>

namespace aspen {

using timer = cpam::timer;

inline auto symmetric_graph_from_static_graph(
    std::tuple<size_t, size_t, uintT*, uintE*>& parsed_graph) {
  using W = empty;
  using inner_graph = symmetric_graph<W>;
  using outer_graph = traversable_graph<inner_graph>;
  timer build_t;
  build_t.start();
  auto G = outer_graph(parsed_graph);
  build_t.stop();
  build_t.reportTotal("Aspen: build time");

  return G;
}

inline auto symmetric_graph_from_static_compressed_graph(char* parsed_graph,
                                                         size_t n_parts = 20) {
  using W = empty;
  using inner_graph = symmetric_graph<W>;
  using outer_graph = traversable_graph<inner_graph>;

  char* s = parsed_graph;
  long* sizes = (long*)s;
  size_t n = sizes[0];
  cout << "sz[0] = " << sizes[0] << " sz[1] = " << sizes[1]
       << " sz[2] = " << sizes[2] << endl;
  uintT* offsets = (uintT*)(s + 3 * sizeof(long));
  long skip = 3 * sizeof(long) + (n + 1) * sizeof(uintT);
  uintE* Degrees = (uintE*)(s + skip);
  skip += n * sizeof(uintE);
  uchar* edges = (uchar*)(s + skip);

  std::cout << "Building compressed graph, nparts = " << n_parts << std::endl;
  size_t vtxs_per_batch = n / n_parts;
  size_t vertices_finished = 0;
  size_t edges_finished = 0;

  timer build_t;
  build_t.start();

  outer_graph G;

  while (vertices_finished < n) {
    size_t bs = std::min(n - vertices_finished, vtxs_per_batch);
    auto offs = parlay::sequence<uintT>::uninitialized(bs);

    parlay::parallel_for(0, bs,
                         [&](size_t i) {
                           size_t degree = Degrees[vertices_finished + i];
                           uintE v = vertices_finished + i;
                           size_t actual_degree = 0;
                           if (degree > 0) {
                             uchar* start = edges + offsets[v];
                             auto it = bytepd_amortized::simple_iter(start,
                                                                     degree, v);
                             uintE ngh = it.cur();
                             if (ngh != v) {
                               actual_degree++;
                             }
                             while (it.has_next()) {
                               uintE nxt = it.next();
                               if (nxt != ngh || nxt != v) {
                                 actual_degree++;
                               }
                               ngh = nxt;
                             }
                           }
                           offs[i] = actual_degree;
                         },
                         1);
    size_t n_edges = parlay::scan_inplace(parlay::make_slice(offs));
    cout << "Next batch num edges = " << n_edges << endl;
    edges_finished += n_edges;

    using ngh_and_weight = typename outer_graph::ngh_and_weight;
    auto batch_edges = parlay::sequence<ngh_and_weight>::uninitialized(n_edges);
    parlay::parallel_for(
        0, bs,
        [&](size_t i) {
          size_t off = offs[i];
          size_t original_degree = Degrees[vertices_finished + i];
          size_t degree = ((i == (bs - 1)) ? (n_edges) : offs[i + 1]) - off;
          if (degree > 0) {
            uintE v = vertices_finished + i;
            uchar* start = edges + offsets[v];
            auto it = bytepd_amortized::simple_iter(start, original_degree, v);
            uintE ngh = it.cur();
            size_t k = 0;
            if (ngh != v) {
              batch_edges[off] = ngh_and_weight(ngh, W());
              k++;
            }
            while (it.has_next()) {
              uintE nxt = it.next();
              if (nxt != ngh || nxt != v) {
                batch_edges[off + k++] = ngh_and_weight(nxt, W());
              }
              ngh = nxt;
            }
          }
        },
        1);

    G.insert_vertex_block(bs, n_edges, offs.begin(), batch_edges.begin(),
                          vertices_finished);

    vertices_finished += bs;
  }
  std::cout << "Finished construction" << std::endl;
  std::cout << "G.n = " << G.num_vertices() << " G.m = " << G.num_edges()
            << std::endl;

  return G;
}

inline auto symmetric_graph_from_static_compressed_graph_edges(
    char* parsed_graph, size_t n_parts = 20) {
  using W = empty;
  using inner_graph = symmetric_graph<W>;
  using outer_graph = traversable_graph<inner_graph>;
  using ngh_and_weight = typename outer_graph::ngh_and_weight;

  char* s = parsed_graph;
  long* sizes = (long*)s;
  size_t n = sizes[0];
  size_t m = sizes[1];
  cout << "sz[0] = " << sizes[0] << " sz[1] = " << sizes[1]
       << " sz[2] = " << sizes[2] << endl;
  uintT* offsets = (uintT*)(s + 3 * sizeof(long));
  long skip = 3 * sizeof(long) + (n + 1) * sizeof(uintT);
  uintE* Degrees = (uintE*)(s + skip);
  skip += n * sizeof(uintE);
  uchar* edges = (uchar*)(s + skip);

  std::cout << "Building compressed graph, nparts = " << n_parts << std::endl;
  size_t edges_per_batch = m / n_parts;
  size_t vertices_finished = 0;
  size_t edges_finished = 0;

  timer build_t;
  build_t.start();

  outer_graph G;
  size_t edges_per_batch_size = edges_per_batch + n;
  auto batch_edges =
      parlay::sequence<ngh_and_weight>::uninitialized(edges_per_batch_size);

  auto scan_degs = parlay::sequence<size_t>(n);
  parlay::parallel_for(0, n, [&](size_t i) { scan_degs[i] = Degrees[i]; });
  size_t total_edges = parlay::scan_inplace(parlay::make_slice(scan_degs));
  std::cout << "m = " << m << " total_edges = " << total_edges << std::endl;

  while (vertices_finished < n) {
    // Compute next batch of vertices.

    size_t next_edges = edges_finished + edges_per_batch;
    size_t vtx_off = parlay::internal::binary_search(
        parlay::make_slice(scan_degs), next_edges, std::less<size_t>());

    size_t bs = vtx_off - vertices_finished;
    auto offs = parlay::sequence<uintT>::uninitialized(bs);

    parlay::parallel_for(0, bs, [&](size_t i) {
      size_t degree = Degrees[vertices_finished + i];
      uintE v = vertices_finished + i;
      size_t actual_degree = 0;
      if (degree > 0) {
        uchar* start = edges + offsets[v];
        auto it = bytepd_amortized::simple_iter(start, degree, v);
        uintE ngh = it.cur();
        if (ngh != v) {
          actual_degree++;
        }
        while (it.has_next()) {
          uintE nxt = it.next();
          if (nxt != ngh || nxt != v) {
            actual_degree++;
          }
          ngh = nxt;
        }
      }
      offs[i] = actual_degree;
    }, 1);
    size_t n_edges = parlay::scan_inplace(parlay::make_slice(offs));
    cout << "Next batch num edges = " << n_edges << endl;

    parlay::parallel_for(0, bs, [&](size_t i) {
      size_t off = offs[i];
      size_t original_degree = Degrees[vertices_finished + i];
      size_t degree = ((i == (bs - 1)) ? (n_edges) : offs[i + 1]) - off;
      if (degree > 0) {
        uintE v = vertices_finished + i;
        uchar* start = edges + offsets[v];
        auto it = bytepd_amortized::simple_iter(start, original_degree, v);
        uintE ngh = it.cur();
        size_t k = 0;
        if (ngh != v) {
          batch_edges[off] = ngh_and_weight(ngh, W());
          k++;
        }
        while (it.has_next()) {
          uintE nxt = it.next();
          if (nxt != ngh || nxt != v) {
            batch_edges[off + k++] = ngh_and_weight(nxt, W());
          }
          ngh = nxt;
        }
      }
    }, 1);

    G.insert_vertex_block(bs, n_edges, offs.begin(), batch_edges.begin(),
                          vertices_finished);

    vertices_finished += bs;
    edges_finished += n_edges;
  }
  std::cout << "Finished construction" << std::endl;
  std::cout << "G.n = " << G.num_vertices() << " G.m = " << G.num_edges()
            << std::endl;

  return G;
}

}  // namespace aspen

#define run_app(G, APP, rounds)             \
  cpam::timer st;                           \
  double total_time = 0.0;                  \
  for (size_t r = 0; r < rounds; r++) {     \
    total_time += APP(G, P);                \
  }                                         \
  auto time_per_iter = total_time / rounds; \
  std::cout << "# time per iter: " << time_per_iter << "\n";

/* Macro to generate binary for unweighted graph applications that can ingest
 * only
 * symmetric graph inputs */
#define generate_symmetric_aspen_main(APP, mutates)                           \
  int main(int argc, char* argv[]) {                                          \
    std::cout << "In main" << std::endl;                                      \
    cpam::commandLine P(argc, argv, " [-s] <inFile>");                        \
    char* iFile = P.getArgument(0);                                           \
    bool symmetric = P.getOptionValue("-s");                                  \
    bool mmap = P.getOptionValue("-m");                                       \
    if (!symmetric) {                                                         \
      std::cout                                                               \
          << "# The application expects the input graph to be symmetric (-s " \
             "flag)."                                                         \
          << std::endl;                                                       \
      std::cout << "# Please run on a symmetric input." << std::endl;         \
    }                                                                         \
    size_t rounds = P.getOptionLongValue("-rounds", 3);                       \
    timer rt;                                                                 \
    rt.start();                                                               \
    bool compressed = P.getOption("-c");                                      \
    if (compressed) {                                                         \
      auto G =                                                                \
          aspen::parse_unweighted_compressed_symmetric_graph(iFile, mmap);    \
      rt.next("Graph read time");                                             \
      auto AG = aspen::symmetric_graph_from_static_compressed_graph_edges(G); \
      run_app(AG, APP, rounds)                                                \
    } else {                                                                  \
      auto G = aspen::parse_unweighted_symmetric_graph(iFile, mmap);          \
      rt.next("Graph read time");                                             \
      auto AG = aspen::symmetric_graph_from_static_graph(G);                  \
      run_app(AG, APP, rounds)                                                \
    }                                                                         \
  }
