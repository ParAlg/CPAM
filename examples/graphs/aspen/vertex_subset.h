#pragma once

#include "flags.h"
#include "macros.h"
#include "utils.h"

#include <functional>
#include <limits>
#include <optional>

namespace aspen {

template <class _data>
struct vertexSubsetData {
  using Data = _data;
  using S = std::tuple<uintE, Data>;
  using D = std::tuple<bool, Data>;

  // An empty vertex set.
  vertexSubsetData(size_t _n) : n(_n), m(0), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // A vertexSubset from array of vertex indices.
  vertexSubsetData(size_t _n, parlay::sequence<S>&& indices)
      : n(_n), m(indices.size()), s_seq(std::move(indices)), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    s = s_seq.begin();
  }

  // A vertexSubset from boolean array giving number of true values.
  vertexSubsetData(size_t _n, size_t _m, parlay::sequence<D>&& _d)
      : n(_n), m(_m), d_seq(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    d = d_seq.begin();
  }

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData(size_t _n, parlay::sequence<D>&& _d) : n(_n), d_seq(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    d = d_seq.begin();
    auto d_map = parlay::delayed_seq<size_t>(n, [&](size_t i) { return (size_t)std::get<0>(_d[i]); });
    m = parlay::reduce(d_map);
  }

  vertexSubsetData() : n(0), m(0), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  bool out_degrees_set() {
    return (sum_out_degrees != std::numeric_limits<size_t>::max());
  }
  size_t get_out_degrees() {
    return sum_out_degrees;
  }
  void set_out_degrees(size_t _sum_out_degrees) {
    sum_out_degrees = _sum_out_degrees;
  }

  bool is_empty() { return m == 0; }

  // Sparse
  inline uintE& vtx(const uintE& i) { return std::get<0>(s[i]); }
  inline Data& vtxData(const uintE& i) { return std::get<1>(s[i]); }
  inline std::tuple<uintE, Data> vtxAndData(const uintE& i) const {
    return s[i];
  }

  // Dense
  __attribute__((always_inline)) inline bool isIn(const uintE& v) const {
    return std::get<0>(d[v]);
  }
  inline Data& ithData(const uintE& v) { return std::get<1>(d[v]); }

  // Returns (uintE) -> std::optional<std::tuple<vertex, vertex-data>>.
  auto get_fn_repr() const
      -> std::function<std::optional<std::tuple<uintE, Data>>(uintE)> {
    std::function<std::optional<std::tuple<uintE, Data>>(const uintE&)> fn;
    if (isDense) {
      fn = [&](const uintE& v) -> std::optional<std::tuple<uintE, Data>> {
        const auto& dv = d[v];
        if (std::get<0>(dv)) {
          return std::optional<std::tuple<uintE, Data>>(std::make_tuple(v, std::get<1>(d[v])));
        } else {
          return std::nullopt;
        }
      };
    } else {
      fn = [&](const uintE& i) -> std::optional<std::tuple<uintE, Data>> {
        return std::optional<std::tuple<uintE, Data>>(s[i]);
      };
    }
    return fn;
  }

  size_t size() const { return m; }
  size_t numVertices() const { return n; }

  size_t numRows() const { return n; }
  size_t numNonzeros() const { return m; }

  bool isEmpty() const { return m == 0; }
  bool dense() const { return isDense; }

  void toSparse() {
    if (s.size() == 0 && m > 0) {
      auto f_seq = parlay::delayed_seq<D>(n, [&](size_t i) -> std::tuple<bool, Data> { return d[i]; });
      auto s_seq = pack_index_and_data<uintE, Data>(f_seq, n);
      if (s_seq.size() != m) {
        abort();
      }
      s = s_seq.begin();
    }
    isDense = false;
  }

  // Convert to dense but keep sparse representation if it exists.
  void toDense() {
    if (d_seq.size() == 0) {
      d_seq = parlay::sequence<D>::uninitialized(n);
      d = d_seq.begin();
      par_for(0, n, [&](size_t i) { std::get<0>(d[i]) = false; });
      par_for(0, m, [&](size_t i) {
        d[std::get<0>(s[i])] = std::make_tuple(true, std::get<1>(s[i]));
      });
    }
    isDense = true;
  }

  size_t n, m;
  parlay::sequence<S> s_seq;
  S* s;
  parlay::sequence<D> d_seq;
  D* d;
  bool isDense;
  size_t sum_out_degrees;
};

// Specialized version where data = empty.
template <>
struct vertexSubsetData<empty> {
  using Data = empty;
  using S = uintE;
  using D = bool;



  // Move constructor
  vertexSubsetData<empty>(vertexSubsetData<empty>&& other) noexcept {
    n = other.n;
    m = other.m;
    s_seq = std::move(other.s_seq);
    s = s_seq.begin();
    d_seq = std::move(other.d_seq);
    d = d_seq.begin();
    isDense = other.isDense;
    sum_out_degrees = other.sum_out_degrees;
  }

  // Move assignment
  vertexSubsetData<empty>& operator=(vertexSubsetData<empty>&& other) noexcept {
    if (this != &other) {
      n = other.n;
      m = other.m;
      s_seq = std::move(other.s_seq);
      s = s_seq.begin();
      d_seq = std::move(other.d_seq);
      d = d_seq.begin();
      isDense = other.isDense;
      sum_out_degrees = other.sum_out_degrees;
    }
    return *this;
  }


  vertexSubsetData<empty>() : n(0), m(0), isDense(0),
      sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // An empty vertex set.
  vertexSubsetData<empty>(size_t _n)
      : n(_n), m(0), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {}

  // A vertexSubset with a single vertex.
  vertexSubsetData<empty>(size_t _n, uintE v)
      : n(_n), m(1), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    s_seq = parlay::sequence<uintE>(1);
    s = s_seq.begin();
    s[0] = v;
  }

  // A vertexSubset from array of vertex indices.
  vertexSubsetData<empty>(size_t _n, parlay::sequence<S>&& indices)
      : n(_n), m(indices.size()), s_seq(std::move(indices)), isDense(0), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    s = s_seq.begin();
  }

  // A vertexSubset from boolean array giving number of true values.
  vertexSubsetData<empty>(size_t _n, size_t _m, parlay::sequence<bool>&& _d)
      : n(_n), m(_m), d_seq(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    d = d_seq.begin();
  }

  // A vertexSubset from boolean array giving number of true values. Calculate
  // number of nonzeros and store in m.
  vertexSubsetData<empty>(size_t _n, parlay::sequence<bool>&& _d)
      : n(_n), d_seq(std::move(_d)), isDense(1), sum_out_degrees(std::numeric_limits<size_t>::max()) {
    d = d_seq.begin();
    auto d_map = parlay::delayed_seq<size_t>(n, [&](size_t i) { return d[i]; });
    m = reduce(d_map);
  }

  bool is_empty() { return m == 0; }

  bool out_degrees_set() {
    return (sum_out_degrees != std::numeric_limits<size_t>::max());
  }
  size_t get_out_degrees() {
    return sum_out_degrees;
  }
  void set_out_degrees(size_t _sum_out_degrees) {
    sum_out_degrees = _sum_out_degrees;
  }

  // Sparse
  inline uintE vtx(const uintE& i) const { return s[i]; }

  inline empty vtxData(const uintE& i) const {
    return empty();
  }
  inline std::tuple<uintE, empty> vtxAndData(const uintE& i) const {
    return std::make_tuple(s[i], empty());
  }

  // Dense
  __attribute__((always_inline)) inline bool isIn(const uintE& v) const {
    return d[v];
  }
  inline empty ithData(const uintE& v) const {
    return empty();
  }

  // Returns (uintE) -> std::optional<std::tuple<vertex, vertex-data>>.
  auto get_fn_repr() const
      -> std::function<std::optional<std::tuple<uintE, empty>>(uintE)> {
    std::function<std::optional<std::tuple<uintE, empty>>(const uintE&)> fn;
    if (isDense) {
      fn = [&](const uintE& v) -> std::optional<std::tuple<uintE, empty>> {
        if (d[v]) {
          return std::optional<std::tuple<uintE, empty>>(std::make_tuple(v, empty()));
        } else {
          return std::nullopt;
        }
      };
    } else {
      fn = [&](const uintE& i) -> std::optional<std::tuple<uintE, empty>> {
        return std::optional<std::tuple<uintE, empty>>(
            std::make_tuple(s[i], empty()));
      };
    }
    return fn;
  }

  size_t size() const { return m; }
  size_t numVertices() const { return n; }

  size_t numRows() const { return n; }
  size_t numNonzeros() const { return m; }

  bool isEmpty() const { return m == 0; }
  bool dense() const { return isDense; }

  void toSparse() {
    if (s_seq.size() == 0 && m > 0) {
      auto _d = d;
      auto f_in = parlay::delayed_seq<bool>(n, [&](size_t i) { return _d[i]; });
      s_seq = parlay::pack_index<uintE>(f_in);
      s = s_seq.begin();
      if (s_seq.size() != m) {
        abort();
      }
    }
    isDense = false;
  }

  // Converts to dense but keeps sparse representation if it exists.
  void toDense() {
    if (d_seq.size() == 0) {
      cpam::timer t; t.start();
      d_seq = parlay::sequence<bool>::uninitialized(n);
      d = d_seq.begin();
      parlay::parallel_for(0, n, [&](size_t i) { d[i] = 0; });
      parlay::parallel_for(0, m, [&](size_t i) { d[s[i]] = 1; });
    }
    isDense = true;
  }

  size_t n, m;
  parlay::sequence<S> s_seq;
  S* s;
  parlay::sequence<bool> d_seq;
  D* d;
  bool isDense;
  size_t sum_out_degrees;
};

using vertexSubset = vertexSubsetData<empty>;

/* ======================== Functions on VertexSubsets ====================== */

// Takes a vertexSubsetData (with some non-trivial Data) and applies a map
// function f : (uintE x Data) -> void over each vertex in the vertexSubset, in
// parallel.
template <class F, class VS,
          typename std::enable_if<!std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f, size_t granularity=2048) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    parlay::parallel_for(0, n, [&] (size_t i) {
      if (V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }, granularity);
  } else {
    parlay::parallel_for(0, m, [&] (size_t i) { f(V.vtx(i), V.vtxData(i)); }, granularity);
  }
}

// Takes a vertexSubset (with no extra data per-vertex) and applies a map
// function f : uintE -> void over each vertex in the vertexSubset, in
// parallel.
template <class VS, class F,
          typename std::enable_if<std::is_same<VS, vertexSubset>::value,
                                  int>::type = 0>
inline void vertexMap(VS& V, F f, size_t granularity=2048) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (V.dense()) {
    parlay::parallel_for(0, n, [&] (size_t i) {
      if (V.isIn(i)) {
        f(i);
      }
    }, granularity);
  } else {
    parlay::parallel_for(0, m, [&] (size_t i)
                    { f(V.vtx(i)); }, granularity);
  }
}

template <class F, class Data>
inline vertexSubset vertexFilter_dense(vertexSubsetData<Data>& V, F filter, size_t granularity=2048) {
  size_t n = V.numRows();
  V.toDense();
  auto d_out = parlay::sequence<bool>(n);
  parlay::parallel_for(0, n, [&] (size_t i) { d_out[i] = 0; }, granularity);
  parlay::parallel_for(0, n, [&] (size_t i) {
    if constexpr (std::is_same<Data, empty>::value) {
      if (V.isIn(i)) d_out[i] = filter(i);
    } else {
      if (V.isIn(i)) d_out[i] = filter(i, V.ithData(i));
    }
  }, granularity);
  return vertexSubset(n, std::move(d_out));
}

template <class F, class Data>
inline vertexSubset vertexFilter_sparse(vertexSubsetData<Data>& V, F filter, size_t granularity=2048) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  auto bits = parlay::sequence<bool>(m);
  V.toSparse();
  parlay::parallel_for(0, m, [&] (size_t i) {
    uintE v = V.vtx(i);
    if constexpr (std::is_same<Data, empty>::value) {
      bits[i] = filter(v);
    } else {
      bits[i] = filter(v, V.vtxData(i));
    }
  }, granularity);
  auto v_imap = parlay::delayed_seq<uintE>(m, [&](size_t i) { return V.vtx(i); });
  auto bits_m = parlay::delayed_seq<bool>(m, [&](size_t i) { return bits[i]; });
  auto out = parlay::pack(v_imap, bits_m);
  size_t out_size = out.size();
  return vertexSubset(n, std::move(out));
}

// Note that this currently strips vertices of their associated data (which is
// the intended use-case in all current uses). Should refactor at some point to
// make keeping/removing the data a choice.
template <class F, class VS>
inline vertexSubset vertexFilter(VS& vs, F filter, flags fl = 0) {
  if (fl == dense_only) {
    return vertexFilter_dense(vs, filter);
  } else if (fl == no_dense) {
    return vertexFilter_sparse(vs, filter);
  }
  if (vs.dense()) {
    return vertexFilter_dense(vs, filter);
  }
  return vertexFilter_sparse(vs, filter);
}

template <class VS>
void add_to_vsubset(VS& vs, uintE* new_verts, uintE num_new_verts) {
  if (vs.dense()) {
    parlay::parallel_for(0, num_new_verts, [&] (size_t i) {
      uintE v = new_verts[i];
      vs.d[v] = true;
    });
  } else {
    vs.s_seq.append(new_verts, new_verts + num_new_verts);
    vs.s = vs.s_seq.begin();
  }
  vs.m += num_new_verts;
}

}  // namespace gbbs
