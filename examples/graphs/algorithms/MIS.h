// This code is part of the project "Theoretically Efficient Parallel Graph
// Algorithms Can Be Fast and Scalable", presented at Symposium on Parallelism
// in Algorithms and Architectures, 2018.
// Copyright (c) 2018 Laxman Dhulipala, Guy Blelloch, and Julian Shun
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all  copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once

#include "aspen/aspen.h"

#include "parlay/random.h"

namespace aspen {
namespace MaximalIndependentSet_rootset {

template <class Graph, class Fl>
inline void verify_mis(Graph& G, Fl& in_mis) {
  using W = typename Graph::weight_type;
  auto d = parlay::sequence<uintE>(G.n, (uintE)0);
  auto map_f = [&](const uintE& src, const uintE& ngh, const W& wgh) {
    if (!d[ngh]) {
      d[ngh] = 1;
    }
  };
  par_for(0, G.n, [&] (size_t i) {
    if (in_mis[i]) {
      G.get_vertex(i).out_neighbors().map(map_f);
    }
  });
  par_for(0, G.n, [&] (size_t i) {
    if (in_mis[i]) {
      assert(!d[i]);
    }
  });
  auto mis_f = [&](size_t i) { return (size_t)in_mis[i]; };
  auto mis_int =
      parlay::delayed_seq<size_t>(G.n, mis_f);
  size_t mis_size = parlay::reduce(mis_int);
  if (parlay::reduce(d) != (G.n - mis_size)) {
    std::cout << "MaximalIndependentSet incorrect"
              << "\n";
    assert(false);
  }
  std::cout << "MaximalIndependentSet Ok"
            << "\n";
}

template <class P, class W>
struct GetNghs {
  P& p;
  GetNghs(P& p) : p(p) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (p[d] > 0) {
      p[d] = 0;
      return true;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    auto p_d = p[d];
    if (p_d > 0 && cpam::utils::atomic_compare_and_swap(&p[d], p_d, 0)) {
      return true;
    }
    return false;
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <class Graph, class VS, class P, class FS>
inline vertexSubset get_nghs(Graph& G, VS& vs, P& p, FS& fs) {
  using W = typename Graph::weight_type;
  if (fs.size() == 0) {
    return G.edgeMap(vs, GetNghs<P, W>(p), -1, no_dense);
  } else {
    return G.edgeMap(vs, GetNghs<P, W>(p), fs, -1, no_dense);
  }
}

inline bool hash_lt(const uintE& src, const uintE& ngh) {
  uint32_t src_h = parlay::hash32(src);
  uint32_t ngh_h = parlay::hash32(ngh);
  return (src_h < ngh_h) || ((src_h == ngh_h) && src < ngh);
};

template <class W>
struct mis_f {
  intE* p;
  uintE* perm;
  mis_f(intE* _p, uintE* _perm) : p(_p), perm(_perm) {}
  inline bool update(const uintE& s, const uintE& d, const W& w) {
    if (perm[s] < perm[d]) {
      p[d]--;
      return p[d] == 0;
    }
    return false;
  }
  inline bool updateAtomic(const uintE& s, const uintE& d, const W& wgh) {
    if (perm[s] < perm[d]) {
      return (cpam::utils::fetch_and_add(&p[d], -1) == 1);
    }
    return false;
  }
  inline bool cond(uintE d) { return (p[d] > 0); }
};

template <class Graph>
inline parlay::sequence<bool> MaximalIndependentSet(Graph& G, bool flatsnap=false) {
  using W = typename Graph::weight_type;
  timer init_t; init_t.start();
  size_t n = G.num_vertices();
  using edge_node = typename Graph::edge_node;
  parlay::sequence<edge_node*> fs;
  if (flatsnap) {
    timer st;
    st.start();
    fs = G.fetch_all_vertices();
    st.next("Snapshot time");
  }

  // compute the priority DAG
  auto priorities = parlay::sequence<intE>::uninitialized(n);
  auto perm = parlay::random_permutation<uintE>(n);

  auto degrees = parlay::sequence<uintE>::uninitialized(n);

  parlay::parallel_for(0, n, [&] (size_t i) {
  //for (size_t i=0; i<n; i++) {
    uintE our_pri = perm[i];
    auto count_f = [&](uintE src, uintE ngh, const W& wgh) {
      uintE ngh_pri = perm[ngh];
      return ngh_pri < our_pri;
    };
    auto vtx = G.get_vertex(i);
    size_t count = 0;
    count = vtx.out_neighbors().count(count_f);
    degrees[i] = vtx.out_degree();
    priorities[i] = count;
  }, 1);
  init_t.stop();
  init_t.reportTotal("init");

  // compute the initial rootset
  auto zero_f = [&](size_t i) { return priorities[i] == 0; };
  auto zero_map =
      parlay::delayed_seq<bool>(n, zero_f);
  auto init = parlay::pack_index<uintE>(zero_map);
  auto roots = vertexSubset(n, std::move(init));

  auto in_mis = parlay::sequence<bool>(n, false);
  size_t finished = 0;
  size_t rounds = 0;
  while (finished != n) {
    assert(roots.size() > 0);
    std::cout << "## round = " << rounds << " size = " << roots.size()
              << " remaining = " << (n - finished) << "\n";

    // set the roots in the MaximalIndependentSet
    vertexMap(roots, [&](uintE v) { in_mis[v] = true; });

    // compute neighbors of roots that are still live using nghMap
    auto removed = get_nghs(G, roots, priorities, fs);
    std::cout << "## removed: " << removed.size() << " many vertices" << std::endl;
    vertexMap(removed, [&](uintE v) { priorities[v] = 0; });

    // compute the new roots: neighbors of removed that have their priorities
    // set to 0 after eliminating all nodes in removed
    intE* pri = priorities.begin();
    timer nr;
    nr.start();
    vertexSubset new_roots;
    if (!flatsnap) {
      new_roots =
        G.edgeMap(removed, mis_f<W>(pri, perm.begin()));
    } else {
      new_roots =
        G.edgeMap(removed, mis_f<W>(pri, perm.begin()), fs);
    }
    nr.stop();
    nr.reportTotal("## new roots time");

    // update finished with roots and removed. update roots.
    finished += roots.size();
    finished += removed.size();

    roots = std::move(new_roots);
    rounds++;
  }
  return in_mis;
}
}  // namespace MaximalIndependentSet_rootset

template <class Graph, class Seq>
inline void verify_MaximalIndependentSet(Graph& G, Seq& mis) {
  using W = typename Graph::weight_type;
  size_t n = G.num_vertices();
  auto ok = parlay::sequence<bool>::from_function(n, [&](size_t i) { return 1; });
  parlay::parallel_for(0, n, [&] (size_t i) {
    auto pred = [&](const uintE& src, const uintE& ngh, const W& wgh) {
      return mis[ngh];
    };
    auto vtx = G.get_vertex(i);
    size_t ct = vtx.out_neighbors().count(pred);
    ok[i] = (mis[i]) ? (ct == 0) : (ct > 0);
  });
  auto ok_f = [&](size_t i) { return ok[i]; };
  auto ok_imap = parlay::delayed_seq<size_t>(n, ok_f);
  size_t n_ok = parlay::reduce(ok_imap);
  if (n_ok == n) {
    std::cout << "valid MaximalIndependentSet"
              << "\n";
  } else {
    std::cout << "invalid MaximalIndependentSet, " << (n - n_ok)
              << " vertices saw bad neighborhoods"
              << "\n";
  }
}

}  // namespace aspen
