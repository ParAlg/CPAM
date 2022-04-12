#pragma once

#include "aspen/aspen.h"

namespace aspen {

template <class W>
struct BFS_F {
  uintE* Parents;
  BFS_F(uintE* _Parents) : Parents(_Parents) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) const {
    Parents[d] = s;
    return 1;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& w) const {
    return (cpam::utils::atomic_compare_and_swap(&Parents[d], UINT_E_MAX, s));
  }
  inline bool cond(const uintE& d) const { return (Parents[d] == UINT_E_MAX); }
};

template <class Graph>
inline parlay::sequence<uintE> BFS(Graph& G, uintE src, bool
    flatsnap=false) {
  using W = typename Graph::weight_type;
  size_t n = G.num_vertices();
  using edge_node = typename Graph::edge_node;

  parlay::sequence<edge_node*> fs;
  if (flatsnap) {
    timer st;
    st.start();
    fs = G.fetch_all_vertices();
    st.next("Snapshot time");
  }

  auto Parents = parlay::sequence<uintE>(n, UINT_E_MAX);
  Parents[src] = src;

  vertexSubset Frontier(n, src);
  size_t reachable = 0;
  while (!Frontier.isEmpty()) {
    std::cout << Frontier.size() << "\n";
    reachable += Frontier.size();
    timer emt;
    emt.start();
    vertexSubset output;
    if (!flatsnap) {
      output = G.edgeMap(Frontier, BFS_F<W>(Parents.begin()), -1,
                                    sparse_blocked | dense_parallel);
    } else {
      output = G.edgeMap(Frontier, BFS_F<W>(Parents.begin()), fs, -1,
                                    sparse_blocked | dense_parallel);
    }
    emt.stop();
    emt.reportTotal("edge map time");
    Frontier = std::move(output);
  }
  std::cout << "Reachable: " << reachable << "\n";
  return Parents;
}

}  // namespace aspen
