#pragma once

#include "aspen/aspen.h"
#include <vector>

namespace aspen {

typedef double fType;

  template <class S, class V, class W>
  struct BC_F {
    S& Scores;
    V& Visited;
    BC_F(S& _Scores, V& _Visited) : Scores(_Scores), Visited(_Visited) {}
    inline bool update(const uintE& s, const uintE& d, const W& w) {
      fType oldV = Scores[d];
      Scores[d] += Scores[s];
      return oldV == 0.0;
    }
    inline bool updateAtomic (const uintE& s, const uintE& d, const W& w) {
      fType to_add = Scores[s];
      fType o_val = cpam::utils::fetch_and_add(&Scores[d],to_add);
      return (o_val == 0);
    }
    inline bool cond (uintE d) { return Visited[d] == 0; }
  };

  template <class W, class S, class V>
  auto make_bc_f(S& scores, V& visited) {
    return BC_F<S, V, W>(scores, visited);
  }

  // Vertex map function to mark visited vertexSubset
  template <class V>
  struct BC_Vertex_F {
    V& Visited;
    BC_Vertex_F(V& _Visited) : Visited(_Visited) {}
    inline bool operator() (uintE i) {
      Visited[i] = 1;
      return 1;
    }
  };

  template <class V>
  auto make_bc_vertex_f(V& visited) {
    return BC_Vertex_F<V>(visited);
  }

  // Vertex map function (used on backwards phase) to mark visited vertexSubset
  // and add to Dependencies score
  template <class V, class D>
  struct BC_Back_Vertex_F {
    V& Visited;
    D& Dependencies;
    D& NumPaths;
    BC_Back_Vertex_F(V& _Visited, D& _Dependencies, D& _NumPaths) :
      Visited(_Visited), Dependencies(_Dependencies), NumPaths(_NumPaths) {}
    inline bool operator() (uintE i) {
      Visited[i] = 1;
      Dependencies[i] += NumPaths[i];
      return 1;
    }
  };

  template <class V, class D>
  auto make_bc_back_vertex_f(V& visited, D& dependencies, D& num_paths) {
    return BC_Back_Vertex_F<V, D>(visited, dependencies, num_paths);
  }

  template <class Graph>
  auto BC(Graph& G, const uintE& start, bool use_flatsnap=true) {
    using W = typename Graph::weight_type;
    size_t n = G.num_vertices();

    using edge_node = typename Graph::edge_node;
    parlay::sequence<edge_node*> fs;
    if (use_flatsnap) {
      timer t; t.start();
      fs = G.fetch_all_vertices();
      t.next("Snapshot time");
    }

    auto NumPaths = parlay::sequence<fType>(n, 0.0);
    auto Visited = parlay::sequence<bool>(n, false);
    Visited[start] = 1; NumPaths[start] = 1.0;

    vertexSubset Frontier(n,start);

    vector<vertexSubset> Levels;

    flags fl = fine_parallel;

    timer sparse_t, dense_t, other_t;
    long round = 0;
    while (!Frontier.is_empty()) {
      round++;
      // std::cout << Frontier.size() << std::endl;
      vertexSubset output = G.edgeMap(Frontier,
          make_bc_f<W>(NumPaths,Visited), fs, -1, fl);
      vertexMap(output, make_bc_vertex_f(Visited)); // mark visited

      Frontier.toSparse();
      Levels.push_back(std::move(Frontier)); // save frontier
      Frontier = std::move(output);
    }
    Levels.push_back(std::move(Frontier));

    auto Dependencies = parlay::sequence<fType>(n, 0.0);

    // Invert numpaths
    parlay::parallel_for(0, n, [&] (size_t i) { NumPaths[i] = 1/NumPaths[i]; });

    parlay::parallel_for(0, n, [&] (size_t i) { Visited[i] = 0; });
    Frontier = std::move(Levels[round-1]);
    vertexMap(Frontier,make_bc_back_vertex_f(Visited,Dependencies,NumPaths));

    std::cout << "Starting reverse. " << std::endl;

    for(long r=round-2;r>=0;r--) {
      G.edgeMap(Frontier, make_bc_f<W>(Dependencies,Visited), fs,
          -1, no_output | fl);
      Frontier = std::move(Levels[r]);
      vertexMap(Frontier,make_bc_back_vertex_f(Visited,Dependencies,NumPaths));
    }

    // Update dependencies scores
    parlay::parallel_for(0, n, [&] (size_t i) {
      Dependencies[i] = (Dependencies[i]-NumPaths[i])/NumPaths[i];
    });
    if (false) {
      sparse_t.report(sparse_t.get_total(), "sparse time");
      dense_t.report(dense_t.get_total(), "dense time");
      other_t.report(other_t.get_total(), "other time");
    }
    return Dependencies;
  }

}  // namespace aspen
