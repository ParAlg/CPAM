#pragma once

#include "flags.h"
#include "vertex_subset.h"

namespace aspen {

// Standard version of edgeMapDense.
template <typename data, typename std::enable_if<
                             std::is_same<data, empty>::value, int>::type = 0>
inline auto get_emdense_gen(bool* next) {
  return [next](uintE ngh, bool m = false) __attribute__((always_inline)) {
    if (m) next[ngh] = 1;
  };
}

template <typename data, typename std::enable_if<
                             !std::is_same<data, empty>::value, int>::type = 0>
inline auto get_emdense_gen(std::tuple<bool, data>* next) {
  return [next](uintE ngh, std::optional<data> m = std::nullopt)
      __attribute__((always_inline)) {
    if (m.has_value()) next[ngh] = std::make_tuple(1, *m);
  };
}

// Standard version of edgeMapSparse.
template <typename data, typename std::enable_if<
                             std::is_same<data, empty>::value, int>::type = 0>
inline auto get_emsparse_gen_full(uintE* outEdges) {
  return [outEdges](uintE ngh, uintT offset, bool m)
      __attribute__((always_inline)) {
    if (m) {
      outEdges[offset] = ngh;
    } else {
      outEdges[offset] = UINT_E_MAX;
    }
  };
}

template <typename data, typename std::enable_if<
                             !std::is_same<data, empty>::value, int>::type = 0>
inline auto get_emsparse_gen_full(std::tuple<uintE, data>* outEdges) {
  return [outEdges](uintE ngh, uintT offset,
                    std::optional<data> m = std::nullopt)
      __attribute__((always_inline)) {
    if (m.has_value()) {
      outEdges[offset] = std::make_tuple(ngh, *m);
    } else {
      std::get<0>(outEdges[offset]) = UINT_E_MAX;
    }
  };
}

// Defines edge_map, vertex_map, adds high level traversal primitives over the
// underlying graph
template <class graph>
struct traversable_graph : private graph {
  using G = graph;
  using G::G;  // imports graph's constructors

  using vertex_tree = typename G::vertex_tree;
  using vertex_node = typename G::vertex_node;
  using vertex_gc = typename G::vertex_gc;
  using edge_tree = typename G::edge_tree;
  using edge_node = typename G::edge_node;
  using vertex = typename G::vertex;
  using weight_type = typename graph::weight_type;
  using ngh_and_weight = typename graph::ngh_and_weight;
  using Empty = empty;

  // for coercing the underlying graph to an traversable_graph
  traversable_graph(graph&& m) {
    if (this != &m) {
      this->V.root = m.V.root;
      m.V.root = nullptr;
    }
  }
  traversable_graph() {}

  template <class Data, class VS, class F>
  auto edgeMapSparse(VS& vs, parlay::sequence<vertex>& vertices, F& f,
                     const flags& fl) {
    using S = typename vertexSubsetData<Data>::S;
    size_t n = num_vertices();

    if (should_output(fl)) {
      auto offsets = parlay::sequence<edge_id>::from_function(
          vertices.size(), [&](size_t i) {
            return (fl & in_edges) ? vertices[i].in_degree()
                                   : vertices[i].out_degree();
          });

      size_t outEdgeCount = parlay::scan_inplace(parlay::make_slice(offsets));
      auto outEdges = parlay::sequence<S>::uninitialized(outEdgeCount);

      parlay::parallel_for(0, vertices.size(), [&](size_t i) {
        edge_id o = offsets[i];
        auto& vtx = vertices[i];
        auto neighbors =
            (fl & in_edges) ? vtx.out_neighbors() : vtx.in_neighbors();
        S* out_edges = outEdges.begin();
        auto g = get_emsparse_gen_full<Data>(out_edges);

        auto map_f = [&](vertex_id v, vertex_id u, auto wgh, size_t i) {
          auto m = f.updateAtomic(v, u, wgh);
          g(u, o + i, m);
        };
        neighbors.map_index(map_f);
      }, 1);
      auto filtered = parlay::filter(
          outEdges, [&](const auto& e) { return e != UINT_E_MAX; });
      return vertexSubsetData<Data>(n, std::move(filtered));
    } else {
      parlay::parallel_for(0, vertices.size(), [&](size_t i) {
        auto& vtx = vertices[i];
        auto neighbors =
            (fl & in_edges) ? vtx.out_neighbors() : vtx.in_neighbors();

        auto map_f = [&](vertex_id v, vertex_id u, auto wgh, size_t i) {
          f.updateAtomic(v, u, wgh);
        };
        neighbors.map_index(map_f);
      }, 1);
      return vertexSubsetData<Data>(n);
    }
  }

  template <class Data, class VS, class F>
  auto edgeMapSparse(VS& vs, parlay::sequence<edge_node*>& flat_snap, F& f,
                     const flags& fl) {
    using S = typename vertexSubsetData<Data>::S;
    size_t n = num_vertices();

    if (should_output(fl)) {
      auto offsets = parlay::sequence<edge_id>::from_function(
          vs.size(), [&](size_t i) {
            uintE v = vs.vtx(i);
            auto vtx = vertex(v, flat_snap[v]);
            return (fl & in_edges) ? vtx.in_degree() : vtx.out_degree();
          });
      size_t outEdgeCount = parlay::scan_inplace(parlay::make_slice(offsets));
      auto outEdges = parlay::sequence<S>::uninitialized(outEdgeCount);
      parlay::parallel_for(0, vs.size(), [&](size_t i) {
        uintE v = vs.vtx(i);
        auto vtx = vertex(v, flat_snap[v]);
        edge_id o = offsets[i];
        auto neighbors =
            (fl & in_edges) ? vtx.out_neighbors() : vtx.in_neighbors();
        S* out_edges = outEdges.begin();
        auto g = get_emsparse_gen_full<Data>(out_edges);

        auto map_f = [&](vertex_id v, vertex_id u, auto wgh, size_t i) {
          auto m = f.updateAtomic(v, u, wgh);
          g(u, o + i, m);
        };
        neighbors.map_index(map_f);
      }, 1);
      auto filtered = parlay::filter(
          outEdges, [&](const auto& e) { return e != UINT_E_MAX; });
      return vertexSubsetData<Data>(n, std::move(filtered));
    } else {
      parlay::parallel_for(0, vs.size(), [&](size_t i) {
        uintE v = vs.vtx(i);
        auto vtx = vertex(v, flat_snap[v]);
        auto neighbors =
            (fl & in_edges) ? vtx.out_neighbors() : vtx.in_neighbors();

        auto map_f = [&](vertex_id v, vertex_id u, auto wgh, size_t i) {
          f.updateAtomic(v, u, wgh);
        };
        neighbors.map_index(map_f);
      }, 1);
      return vertexSubsetData<Data>(n);
    }
  }

  template <class F>
  vertexSubset edgeMapDense(vertexSubset& vs, F& f, const flags& fl) {
    size_t n = num_vertices();
    vs.toDense();
    auto prev = vs.d;

    if (should_output(fl)) {
      auto next = parlay::sequence<bool>(n, false);
      auto map_f = [&](const auto& vtx) {
        vertex_id v = vtx.id;
        if (f.cond(v)) {
          auto neighbors = vtx.in_neighbors();
          auto map_cond = [&](const vertex_id& v, const vertex_id& ngh,
                              const auto& wgh) -> bool {
            if (prev[ngh] && f.update(ngh, v, wgh)) {
              next[v] = true;
            }
            return f.cond(v);
          };
          // auto cond = [&] () { return f.cond(v); };
          // neighbors.map_cond(map_cond, cond);
          neighbors.foreach_cond(map_cond);
        }
      };
      map_vertices(map_f);
      return vertexSubset(n, std::move(next));
    } else {
      auto map_f = [&](const auto& vtx) {
        vertex_id v = vtx.id;
        if (f.cond(v)) {
          auto neighbors = vtx.in_neighbors();
          auto map_cond = [&](const vertex_id& v, const vertex_id& ngh,
                              const auto& wgh) -> bool {
            if (prev[ngh]) {
              f.update(ngh, v, wgh);
            }
            return f.cond(v);
          };
          neighbors.foreach_cond(map_cond);
        }
      };
      map_vertices(map_f);
      return vertexSubset(n);
    }
  }

  template <class F>
  vertexSubset edgeMapDense(vertexSubset& vs,
                            parlay::sequence<edge_node*>& flat_snap, F& f,
                            const flags& fl) {
    size_t n = num_vertices();
    vs.toDense();

    auto prev = vs.d;

    if (should_output(fl)) {
      auto next = parlay::sequence<bool>(n, false);
      size_t granularity = (fl & fine_parallel) ? 1 : 1024;
      parlay::parallel_for(0, n, [&] (size_t v) {
        if (f.cond(v)) {
          auto vtx = vertex(v, flat_snap[v]);
          auto neighbors = vtx.in_neighbors();
          auto map_cond = [&](const vertex_id& v, const vertex_id& ngh,
                              const auto& wgh) -> bool {
            if (prev[ngh] && f.update(ngh, v, wgh)) {
              next[v] = true;
            }
            return f.cond(v);
          };
          neighbors.foreach_cond(map_cond);
        }
      }, granularity);
      return vertexSubset(n, std::move(next));
    } else {
      size_t granularity = (fl & fine_parallel) ? 1 : 1024;
      parlay::parallel_for(0, n, [&] (size_t v) {
        if (f.cond(v)) {
          auto vtx = vertex(v, flat_snap[v]);
          auto neighbors = vtx.in_neighbors();
          auto map_cond = [&](const vertex_id& v, const vertex_id& ngh,
                              const auto& wgh) -> bool {
            if (prev[ngh]) {
              f.update(ngh, v, wgh);
            }
            return f.cond(v);
          };
          neighbors.foreach_cond(map_cond);
        }
      }, granularity);
      return vertexSubset(n);
    }
  }

  template <class Data, class VS, class F>
  auto edgeMapData(VS& vs, F& f, long threshold = -1, const flags& fl = 0) {
    size_t n = num_vertices();
    size_t m = num_edges();
    if (threshold == -1) threshold = m / 20;

    if (vs.isDense && vs.size() > (n / 10)) {
      return edgeMapDense<F>(vs, f, fl);
    }

    if (vs.size() == 0) return vertexSubsetData<Data>(n);

    auto vertices = parlay::sequence<vertex>::uninitialized(vs.size());
    size_t out_degrees = 0;
    bool init_vertices = false;
    if (vs.out_degrees_set()) {
      out_degrees = vs.get_out_degrees();
    } else {
      vs.toSparse();
      init_vertices = true;
      parlay::parallel_for(
          0, vs.size(), [&](size_t i) { vertices[i] = get_vertex(vs.vtx(i)); },
          1);

      auto degree_f = [&](size_t i) {
        return (fl & in_edges) ? vertices[i].in_degree()
                               : vertices[i].out_degree();
      };
      auto degree_im = parlay::delayed_seq<size_t>(vertices.size(), degree_f);
      out_degrees = parlay::reduce(degree_im);
      vs.set_out_degrees(out_degrees);
    }
    if (out_degrees == 0) return vertexSubsetData<Data>(n);
    if (vs.size() + out_degrees > static_cast<size_t>(threshold) &&
        !(fl & no_dense)) {
      return edgeMapDense<F>(vs, f, fl);
    }
    if (!init_vertices) {
      parlay::parallel_for(0, vs.size(), [&](size_t i) { vertices[i] = get_vertex(vs.vtx(i)); }, 1);
    }
    return edgeMapSparse<Data, VS, F>(vs, vertices, f, fl);
  }

  template <class VS, class F>
  auto edgeMap(VS& vs, F f, long threshold = -1, const flags& fl = 0) {
    return edgeMapData<Empty, VS, F>(vs, f, threshold, fl);
  }

  template <class Data, class VS, class F>
  auto edgeMapData(VS& vs, F& f, parlay::sequence<edge_node*>& flat_snap,
                   long threshold = -1, const flags& fl = 0) {
    size_t n = num_vertices();
    size_t m = num_edges();
    if (threshold == -1) threshold = m / 20;

    if (vs.isDense && vs.size() > (n / 10)) {
      return edgeMapDense<F>(vs, flat_snap, f, fl);
    }

    if (vs.size() == 0) return vertexSubsetData<Data>(n);

    size_t out_degrees = 0;
    if (vs.out_degrees_set()) {
      out_degrees = vs.get_out_degrees();
    } else {
      vs.toSparse();
      auto degree_f = [&](size_t i) {
        uintE v = vs.vtx(i);
        auto vtx = vertex(v, flat_snap[v]);
        return (fl & in_edges) ? vtx.in_degree() : vtx.out_degree();
      };
      auto degree_im = parlay::delayed_seq<size_t>(vs.size(), degree_f);
      out_degrees = parlay::reduce(degree_im);
      vs.set_out_degrees(out_degrees);
    }
    if (out_degrees == 0) return vertexSubsetData<Data>(n);
    if (vs.size() + out_degrees > static_cast<size_t>(threshold) &&
        !(fl & no_dense)) {
      return edgeMapDense<F>(vs, flat_snap, f, fl);
    }
    return edgeMapSparse<Data, VS, F>(vs, flat_snap, f, fl);
  }

  template <class VS, class F>
  auto edgeMap(VS& vs, F f, parlay::sequence<edge_node*>& flat_snap,
               long threshold = -1, const flags& fl = 0) {
    if (flat_snap.size() == 0) {
      return edgeMapData<Empty, VS, F>(vs, f, threshold, fl);
    } else {
      return edgeMapData<Empty, VS, F>(vs, f, flat_snap, threshold, fl);
    }
  }

  parlay::sequence<edge_node*> fetch_all_vertices() {
    size_t n = G::num_vertices();
    auto vtxs = parlay::sequence<edge_node*>(n, nullptr);
    auto map_f = [&](const auto& vtx) {
      const uintE& v = vtx.id;
      vtxs[v] = vtx.edges;
    };
    G::map_vertices(map_f);
    return vtxs;
  }

  template <class Edge>
  traversable_graph insert_edges_batch(size_t m, Edge* edges, bool sorted = false,
                          bool remove_dups = false,
                          size_t nn = std::numeric_limits<size_t>::max(),
                          bool run_seq = false) {
    return traversable_graph(G::insert_edges_batch_2(m, edges, sorted, remove_dups, nn, run_seq));
  }

  template <class Edge>
  traversable_graph delete_edges_batch(size_t m, Edge* edges, bool sorted = false,
                          bool remove_dups = false,
                          size_t nn = std::numeric_limits<size_t>::max(),
                          bool run_seq = false) {
    // G::delete_edges_batch_1(m, edges, sorted, remove_dups, nn, run_seq);
    return traversable_graph(G::delete_edges_batch_2(m, edges, sorted, remove_dups, nn, run_seq));
  }


  template <class Edge>
  void insert_edges_batch_inplace(size_t m, Edge* edges, bool sorted = false,
                          bool remove_dups = false,
                          size_t nn = std::numeric_limits<size_t>::max(),
                          bool run_seq = false) {
    G::insert_edges_batch(m, edges, sorted, remove_dups, nn, run_seq);
  }

  template <class Edge>
  void delete_edges_batch_inplace(size_t m, Edge* edges, bool sorted = false,
                          bool remove_dups = false,
                          size_t nn = std::numeric_limits<size_t>::max(),
                          bool run_seq = false) {
    G::delete_edges_batch(m, edges, sorted, remove_dups, nn, run_seq);
  }

  // Only used when building compressed (large) graphs
  void insert_vertex_block(size_t vn, size_t vm, uintT* offsets, ngh_and_weight* edges, size_t vertex_offset) {
    G::insert_vertex_block(vn, vm, offsets, edges, vertex_offset);
  }

  void set_num_vertices(size_t _n) {
    this->n = _n;
  }

  using G::get_tree_sizes;
  using G::print_stats;
  using G::num_vertices;
  using G::num_edges;
  using G::get_vertex;
  using G::get_vertices;  // unused?
  using G::map_vertices;
  using G::get_root;
  using G::clear_root;
  using G::set_root;
  using G::ref_cnt;
  using G::iterate_seq;
};


}  // namespace aspen
