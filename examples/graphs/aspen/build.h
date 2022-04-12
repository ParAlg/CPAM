#pragma once

#include "macros.h"
#include "utils.h"

// #define USE_PAM 1

namespace aspen {

// Builds an aspen graph from the given static graph.
template <class weight>
auto graph_to_edges(std::tuple<size_t, size_t, uintT*, uintE*>& parsed_graph,
    size_t num_batches = 1) {
  size_t n = std::get<0>(parsed_graph);

  auto offsets = std::get<2>(parsed_graph);
  auto E = std::get<3>(parsed_graph);
  auto degs = parlay::sequence<size_t>::from_function(
      n, [&](size_t i) { return offsets[i+1] - offsets[i]; });
  size_t sum_degs = parlay::scan_inplace(parlay::make_slice(degs));
  assert(sum_degs == std::get<1>(parsed_graph));

  #ifdef USE_PAM
  using ngh_and_weight = std::pair<vertex_id, weight>;
  #else
  using ngh_and_weight = std::tuple<vertex_id, weight>;
  #endif
  using edge = std::pair<vertex_id, ngh_and_weight>;


  auto edges = parlay::sequence<edge>::uninitialized(sum_degs);

  parlay::parallel_for(0, n, [&](size_t i) {
    size_t k = degs[i];
    auto deg = offsets[i+1] - offsets[i];
    size_t offset = offsets[i];
    parlay::parallel_for(0, deg, [&] (size_t j) {
      edges[k + j] = std::make_pair(i, ngh_and_weight(E[offset + j], empty()));
    });
  }, 1);
  return edges;
}

}  // namespace aspen
