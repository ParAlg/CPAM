#pragma once

#include <cpam/cpam.h>
#include <pam/get_time.h>
#include <pam/pam.h>
#include <pam/parse_command_line.h>

#include <mutex>

#include "../graphs/aspen/aspen.h"
#include "types.h"
#include "util/NSGDist.h"
#include "util/counter.h"
#include "util/table.h"
#include "flock_private_structs/hash_block/set.h"

bool stats = false;

template <typename T>
struct knn_index {
  // Pointers to the coordinates for each point.
  parlay::sequence<Tvec_point<T>*>& v;

  size_t maxDeg;
  size_t beamSize;
  double alpha;
  size_t d;
  Distance* D;
  atomic_sum_counter<size_t> total_visited;

  struct empty_weight {};
  using Graph = aspen::symmetric_graph<empty_weight>;
  using vertex = typename Graph::vertex;
  using edge_tree = typename Graph::edge_tree;
  using edge_node = typename Graph::edge_node;
  using vertex_tree = typename Graph::vertex_tree;
  using vertex_node = typename Graph::vertex_node;
  using LockGuard = std::lock_guard<std::mutex>;

  aspen::versioned_graph<Graph> VG;

  using tvec_point = Tvec_point<T>;
  using fvec_point = Tvec_point<float>;

  tvec_point* medoid;

  using pid = std::pair<node_id, float>;
  using slice_tvec = decltype(make_slice(parlay::sequence<tvec_point*>()));

  using FS = flock_set<node_id, parlay::slice<T*, T*>>;
  FS data_store;

  knn_index(parlay::sequence<Tvec_point<T>*>& v, size_t maxDeg, size_t beamSize,
            double Alpha, size_t dim, Distance* DD)
      : v(v), maxDeg(maxDeg), beamSize(beamSize), alpha(Alpha), d(dim), D(DD) {
    std::cout << "Initialized knn_index with maxDeg = " << maxDeg
              << " beamSize = " << beamSize << " alpha = " << alpha
              << " dim = " << dim << " distance function: " << D->id()
              << std::endl;

    // FS data_store = flock_set<node_id, parlay::slice<T*, T*>>(v.size());
    timer table_t;
    table_t.start();
    data_store.initialize(v.size()*2);
    parlay::parallel_for(0, v.size(), [&] (size_t i){data_store.insert((node_id) i, v[i]->coordinates);});
    table_t.next("Insert time");
    total_visited.reset();
  }

  void print_graph_status() {
    auto S = VG.acquire_version();
    Graph G = S.graph;
    std::cout << "G.num_vertices = " << G.num_vertices()
              << " num_edges = " << G.num_edges() << std::endl;
    VG.release_version(std::move(S));
  }

  void build_index(parlay::sequence<node_id> inserts) {
    // Find the medoid, which each beamSearch will begin from.
    std::cout << "Distance calls: " << distance_calls.get_value() << std::endl;
    std::cout << "Total visited: " << total_visited.get_value() << std::endl;
    Graph Initial_Graph;
    node_id medoid_id = find_approx_medoid();
    Initial_Graph.insert_vertex_inplace(medoid_id, nullptr);
    batch_insert_inplace(Initial_Graph, inserts, 2, .02, false);
    VG = aspen::versioned_graph<Graph>(std::move(Initial_Graph));
    auto S = VG.acquire_version();
    std::cout << "G.num_vertices = " << S.graph.num_vertices()
              << " num_edges = " << S.graph.num_edges()
              << " timestamp = " << S.timestamp
              << " ref_cnt = " << S.get_ref_cnt()
              << " graph_ref_cnt = " << S.get_graph_ref_cnt() << std::endl;
    VG.release_version(std::move(S));
#ifdef STATS
    std::cout << "Total vertices visited: " << total_visited.get_value()
              << std::endl;
    std::cout << "Total distance calls: " << distance_calls.get_value()
              << std::endl;
#endif
  }

  void insert(node_id p) {
    parlay::sequence<node_id> inserts = {p};
    insert(inserts);
  }

  void insert(parlay::sequence<node_id> inserts, bool random_order = false) {
    auto S = VG.acquire_version();
    std::cout << "Acquired version with timestamp " << S.timestamp
              << " address = " << ((size_t)S.get_root())
              << " ref_cnt = " << S.get_ref_cnt() << std::endl;

    parlay::sequence<node_id> init_inserts({inserts[0]});
    parlay::sequence<node_id> remaining_inserts(inserts.size() - 1);

    parlay::parallel_for(0, inserts.size(), [&](size_t i) {
      if (i == 0) {
        init_inserts[i] = inserts[i];
      } else {
        remaining_inserts[i - 1] = inserts[i];
      }
    });

    std::cout << "Performing functional update: update size = "
              << init_inserts.size() << std::endl;
    Graph new_G = batch_insert_functional(S.graph, inserts);
    std::cout << "Performing inplace update on node: "
              << ((size_t)new_G.get_root()) << std::endl;

    batch_insert_inplace(new_G, remaining_inserts, 2, .02, random_order);
    std::cout << "Posting new version " << std::endl;
    VG.add_version_from_graph(std::move(new_G));
    VG.release_version(std::move(S));
  }

  void query(
      parlay::sequence<Tvec_point<T>*> &q, int k, int beamSizeQ,
      float cut, parlay::sequence<parlay::sequence<node_id>> &neighbors) {
    if ((k + 1) > beamSizeQ) {
      std::cout << "Error: beam search parameter Q = " << beamSizeQ
                << " same size or smaller than k = " << k << std::endl;
      abort();
    }
    auto S = VG.acquire_version();
    parlay::parallel_for(0, q.size(), [&](size_t i) {
      auto pairs =
          beam_search_2(S.graph, q[i]->coordinates.begin(), beamSizeQ, k, cut);
      auto& beamElts = pairs.first;
      parlay::sequence<node_id> single_neighbors(k);
      // Ignoring reporting the point itself for now.
      for (int j = 0; j < k; j++) {
        single_neighbors[j] = beamElts[j].first;
      }
      neighbors[i] = single_neighbors;
    });
    VG.release_version(std::move(S));
  }

  // TODO have multiple queries use the same version
  parlay::sequence<node_id> query(T* query_coords, int k, int beamSizeQ,
                                  float cut) {
    if ((k + 1) > beamSizeQ) {
      std::cout << "Error: beam search parameter Q = " << beamSizeQ
                << " same size or smaller than k = " << k << std::endl;
      abort();
    }
    auto S = VG.acquire_version();
    auto pairs = beam_search(S.graph, query_coords, beamSizeQ);
    auto& beamElts = pairs.first;
    parlay::sequence<node_id> neighbors(k);
    // Ignoring reporting the point itself for now.
    for (int j = 0; j < k; j++) {
      neighbors[j] = beamElts[j].first;
    }
    VG.release_version(std::move(S));
    return std::move(neighbors);
  }

  void lazy_delete(parlay::sequence<node_id> deletes) {
    {
      LockGuard guard(delete_lock);
      for (node_id p : deletes) {
        if (p < 0 || p > (node_id)v.size()) {
          std::cout << "ERROR: invalid point " << p << " given to lazy_delete"
                    << std::endl;
          abort();
        }
        if (p != medoid->id)
          delete_set.insert(p);

        else
          std::cout << "Deleting medoid not permitted; continuing" << std::endl;
      }
    }
  }

  void lazy_delete(node_id p) {
    if (p < 0 || p > (node_id)v.size()) {
      std::cout << "ERROR: invalid point " << p << " given to lazy_delete"
                << std::endl;
      abort();
    }
    if (p == (node_id)medoid->id) {
      std::cout << "Deleting medoid not permitted; continuing" << std::endl;
      return;
    }
    {
      LockGuard guard(delete_lock);
      delete_set.insert(p);
    }
  }

  void start_delete_epoch() {
    // freeze the delete set and start a new one before consolidation
    if (!epoch_running) {
      {
        LockGuard guard(delete_lock);
        delete_set.swap(old_delete_set);
      }
      epoch_running = true;
    } else {
      std::cout
          << "ERROR: cannot start new epoch while previous epoch is running"
          << std::endl;
      abort();
    }
  }

  void end_delete_epoch() {
    if (epoch_running) {
      parlay::sequence<node_id> delete_vec;
      for (auto d : old_delete_set) delete_vec.push_back(d);

      auto W = VG.acquire_version();
      Graph new_G = W.graph.delete_vertices_batch_functional(
          delete_vec.size(), delete_vec.begin());
      check_deletes_correct(new_G);
      VG.add_version_from_graph(new_G);
      VG.release_version(std::move(W));

      old_delete_set.clear();
      epoch_running = false;
    } else {
      std::cout << "ERROR: cannot end epoch while epoch is not running"
                << std::endl;
      abort();
    }
  }

  void consolidate_deletes(parlay::sequence<node_id> to_consolidate) {
    if (epoch_running) {
      auto S = VG.acquire_version();
      // Graph new_G = consolidate_deletes_with_pruning(S.graph,
      // to_consolidate);
      Graph new_G = consolidate_deletes_with_pruning(S.graph, to_consolidate);
      VG.add_version_from_graph(new_G);
      VG.release_version(std::move(S));
    } else {
      std::cout << "ERROR: cannot consolidate when delete epoch not initialized"
                << std::endl;
      abort();
    }
  }

  node_id get_medoid() { return medoid->id; }

  void print_graph_stats(){
    auto S = VG.acquire_version();
    std::string graph = "g";
    std::string mode = "m";
    S.graph.get_tree_sizes(graph, mode);
    VG.release_version(std::move(S));
  }

 private:
  std::set<node_id> delete_set;
  std::set<node_id> old_delete_set;
  std::mutex delete_lock;  // lock for delete_set which can only be updated
                           // sequentially
  bool epoch_running = false;
  // p_coords: query vector coordinates
  // v: database of vectors
  // medoid: "root" of the proximity graph
  // beamSize: (similar to ef)
  // d: dimensionality of the indexed vectors



  Graph consolidate_deletes_simple(Graph G) {
    auto consolidated_vertices =
        parlay::sequence<std::tuple<node_id, edge_node*>>(v.size());
    auto needs_consolidate = parlay::sequence<bool>(v.size(), false);

    parlay::parallel_for(0, v.size(), [&](size_t i) {
      if (old_delete_set.find(i) == old_delete_set.end()) {
        auto current_vtx = G.get_vertex(i);
        parlay::sequence<node_id> candidates;
        auto g = [&](node_id a) {
          return old_delete_set.find(a) == old_delete_set.end();
        };
        auto f = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v)) candidates.push_back(v);
          return true;
        };
        current_vtx.out_neighbors().foreach_cond(f);
        auto begin = (std::tuple<node_id, empty_weight>*)candidates.begin();
        auto tree = edge_tree(begin, begin + candidates.size());
        consolidated_vertices[i] = {i, tree.root};
        needs_consolidate[i] = true;
        tree.root = nullptr;
      }
    });
    auto filtered_vertices =
        parlay::pack(consolidated_vertices, needs_consolidate);
    Graph new_G = G.insert_vertices_batch_functional(filtered_vertices.size(),
                                                     filtered_vertices.begin());
    return new_G;
  }

  Graph consolidate_deletes_with_pruning(
      Graph G, parlay::sequence<node_id>& to_consolidate) {
    auto consolidated_vertices =
        parlay::sequence<std::tuple<node_id, edge_node*>>(
            to_consolidate.size());
    auto needs_consolidate =
        parlay::sequence<bool>(to_consolidate.size(), false);

    parlay::parallel_for(0, to_consolidate.size(), [&](size_t i) {
      node_id index = to_consolidate[i];
      if (old_delete_set.find(index) == old_delete_set.end()) {
        // remove deleted vertices and add their out_nbh
        bool change = false;
        auto current_vtx = G.get_vertex(index);
        parlay::sequence<node_id> candidates;
        auto g = [&](node_id a) {
          return old_delete_set.find(a) == old_delete_set.end();
        };

        auto f = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v)) candidates.push_back(v);
          return true;
        };

        auto h = [&](node_id u, node_id v, empty_weight wgh) {
          if (g(v))
            candidates.push_back(v);
          else {
            change = true;
            auto vtx = G.get_vertex(v);
            vtx.out_neighbors().foreach_cond(f);
          }
          return true;
        };
        current_vtx.out_neighbors().foreach_cond(h);
        if (change) {
          needs_consolidate[i] = true;
          if (candidates.size() <= maxDeg) {
            auto begin = (std::tuple<node_id, empty_weight>*)candidates.begin();
            auto tree = edge_tree(begin, begin + candidates.size());
            consolidated_vertices[i] = {index, tree.root};
            tree.root = nullptr;
          } else {
            parlay::sequence<node_id> new_out_2(
                maxDeg, std::numeric_limits<node_id>::max());
            auto output_slice = parlay::make_slice(new_out_2.begin(),
                                                   new_out_2.begin() + maxDeg);
            robustPrune(G, v[index], index, candidates, alpha, output_slice,
                        false);
            size_t deg = size_of(output_slice);
            auto begin = (std::tuple<node_id, empty_weight>*)new_out_2.begin();
            auto tree = edge_tree(begin, begin + deg);
            consolidated_vertices[i] = {index, tree.root};
            tree.root = nullptr;
          }
        }
      }
    });
    auto filtered_vertices =
        parlay::pack(consolidated_vertices, needs_consolidate);
    Graph new_G = G.insert_vertices_batch_functional(filtered_vertices.size(),
                                                     filtered_vertices.begin());
    return new_G;
  }

  void check_deletes_correct(Graph& G) {
    auto map = [&](vertex current_vtx) {
      auto g = [&](node_id a) {
        return (old_delete_set.find(a) != old_delete_set.end());
      };
      auto f = [&](node_id u, node_id v, empty_weight wgh) {
        if (g(v)) {
          std::cout << "ERROR: vertex " << u << " has deleted neighbor " << v
                    << std::endl;
        }
        return true;
      };
      current_vtx.out_neighbors().foreach_cond(f);
    };
    G.map_vertices(map);
  }

   // updated version by Guy
  std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search_2(
      Graph& G, T* p_coords, int beamSize, int k = 0, float cut = 1.14) {
    // initialize data structures

    std::vector<pid> visited;
    parlay::sequence<pid> frontier;
    auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first);
    };
    auto make_pid = [&](node_id q) -> std::pair<node_id, double> {
      auto dist = D->distance(data_store.find(q).begin(), p_coords, d);
      return std::pair{q, dist};
    };

    int bits = std::ceil(std::log2(beamSize * beamSize)) - 2;
    parlay::sequence<node_id> hash_table(1 << bits,
                                         std::numeric_limits<node_id>::max());

    // the frontier starts with the medoid
    frontier.push_back(make_pid(medoid->id));

    std::vector<pid> unvisited_frontier(beamSize);
    parlay::sequence<pid> new_frontier(beamSize + maxDeg);
    unvisited_frontier[0] = frontier[0];
    int remain = 1;

    // terminate beam search when the entire frontier has been visited
    while (remain > 0) {
      // the next node to visit is the unvisited frontier node that is closest
      // to p
      pid currentPid = unvisited_frontier[0];
      auto current_vtx = G.get_vertex(currentPid.first);

      auto g = [&](node_id a) {
        node_id loc = parlay::hash64_2(a) & ((1 << bits) - 1);
        if (hash_table[loc] == a) return false;
        hash_table[loc] = a;
        return true;
      };

      parlay::sequence<node_id> candidates;
      auto f = [&](node_id u, node_id v, empty_weight wgh) {
        if (g(v)) {
          candidates.push_back(v);
        }
        return true;
      };
      current_vtx.out_neighbors().foreach_cond(f);
      auto pairCandidates =
          parlay::map(candidates, [&](long c) { return make_pid(c); });
      auto sortedCandidates = parlay::sort(pairCandidates, less);
      auto f_iter = std::set_union(
          frontier.begin(), frontier.end(), sortedCandidates.begin(),
          sortedCandidates.end(), new_frontier.begin(), less);
      size_t f_size = std::min<size_t>(beamSize, f_iter - new_frontier.begin());
      if (k > 0 && (int)f_size > k)
        f_size = (std::upper_bound(
                      new_frontier.begin(), new_frontier.begin() + f_size,
                      std::pair{0, cut * new_frontier[k].second}, less) -
                  new_frontier.begin());
      frontier =
          parlay::tabulate(f_size, [&](long i) { return new_frontier[i]; });
      visited.insert(
          std::upper_bound(visited.begin(), visited.end(), currentPid, less),
          currentPid);
      auto uf_iter =
          std::set_difference(frontier.begin(), frontier.end(), visited.begin(),
                              visited.end(), unvisited_frontier.begin(), less);
      remain = uf_iter - unvisited_frontier.begin();
    }
#ifdef STATS
    total_visited.update_value(visited.size());
#endif
    //TODO should we also lock the current delete set and filter those elements out?
    parlay::sequence<pid> to_filter = parlay::to_sequence(visited);
    auto f = [&] (pid a){
      return old_delete_set.find(a.first) == old_delete_set.end();
    };
    auto filtered = parlay::filter(to_filter, f);
    return std::make_pair(frontier, std::move(filtered));
  }

  // updated version by Guy
  std::pair<parlay::sequence<pid>, parlay::sequence<pid>> beam_search(
      Graph& G, T* p_coords, int beamSize, int k = 0, float cut = 1.14) {
    // initialize data structures
    auto vvc = v[0]->coordinates.begin();
    long stride = v[1]->coordinates.begin() - v[0]->coordinates.begin();

    std::vector<pid> visited;
    parlay::sequence<pid> frontier;
    auto less = [&](pid a, pid b) {
      return a.second < b.second || (a.second == b.second && a.first < b.first);
    };
    auto make_pid = [&](node_id q) -> std::pair<node_id, double> {
      auto dist = D->distance(v[q]->coordinates.begin(), p_coords, d);
      return std::pair{q, dist};
    };

    int bits = std::ceil(std::log2(beamSize * beamSize)) - 2;
    parlay::sequence<node_id> hash_table(1 << bits,
                                         std::numeric_limits<node_id>::max());

    // the frontier starts with the medoid
    frontier.push_back(make_pid(medoid->id));

    std::vector<pid> unvisited_frontier(beamSize);
    parlay::sequence<pid> new_frontier(beamSize + maxDeg);
    unvisited_frontier[0] = frontier[0];
    int remain = 1;

    // terminate beam search when the entire frontier has been visited
    while (remain > 0) {
      // the next node to visit is the unvisited frontier node that is closest
      // to p
      pid currentPid = unvisited_frontier[0];
      auto current_vtx = G.get_vertex(currentPid.first);

      auto g = [&](node_id a) {
        node_id loc = parlay::hash64_2(a) & ((1 << bits) - 1);
        if (hash_table[loc] == a) return false;
        hash_table[loc] = a;
        return true;
      };

      parlay::sequence<node_id> candidates;
      auto f = [&](node_id u, node_id v, empty_weight wgh) {
        if (g(v)) {
          candidates.push_back(v);
        }
        return true;
      };
      current_vtx.out_neighbors().foreach_cond(f);
      auto pairCandidates =
          parlay::map(candidates, [&](long c) { return make_pid(c); });
      auto sortedCandidates = parlay::sort(pairCandidates, less);
      auto f_iter = std::set_union(
          frontier.begin(), frontier.end(), sortedCandidates.begin(),
          sortedCandidates.end(), new_frontier.begin(), less);
      size_t f_size = std::min<size_t>(beamSize, f_iter - new_frontier.begin());
      if (k > 0 && (int)f_size > k)
        f_size = (std::upper_bound(
                      new_frontier.begin(), new_frontier.begin() + f_size,
                      std::pair{0, cut * new_frontier[k].second}, less) -
                  new_frontier.begin());
      frontier =
          parlay::tabulate(f_size, [&](long i) { return new_frontier[i]; });
      visited.insert(
          std::upper_bound(visited.begin(), visited.end(), currentPid, less),
          currentPid);
      auto uf_iter =
          std::set_difference(frontier.begin(), frontier.end(), visited.begin(),
                              visited.end(), unvisited_frontier.begin(), less);
      remain = uf_iter - unvisited_frontier.begin();
    }
#ifdef STATS
    total_visited.update_value(visited.size());
#endif
    //TODO should we also lock the current delete set and filter those elements out?
    parlay::sequence<pid> to_filter = parlay::to_sequence(visited);
    auto f = [&] (pid a){
      return old_delete_set.find(a.first) == old_delete_set.end();
    };
    auto filtered = parlay::filter(to_filter, f);
    return std::make_pair(frontier, std::move(filtered));
  }

    // robustPrune routine as found in DiskANN paper.
  // The new candidate set is added to the supplied array (new_nghs).
  void robustPrune(Graph& GG, T* p_coords, node_id p_id,
                   parlay::sequence<pid> candidates, double alpha,
                   parlay::slice<node_id*, node_id*> new_nghs,
                   bool add = true) {
    auto less = [&](pid a, pid b) { return a.second < b.second; };
    // add out neighbors of p to the candidate set.
    if (add) {
      auto p_vertex = GG.get_vertex(p_id);
      auto map_f = [&](node_id u, node_id neighbor_id, empty_weight weight) {
        candidates.push_back(std::make_pair(
            neighbor_id, D->distance(data_store.find(neighbor_id).begin(),
                                     p_coords, d)));
        return true;
      };
      p_vertex.out_neighbors().foreach_cond(map_f);
    }
    if (candidates.size() <= maxDeg) {
      uint32_t offset = 0;
      parlay::sort_inplace(candidates, less);
      for (size_t i = 0; i < candidates.size(); ++i) {
        node_id ngh = candidates[i].first;
        if (ngh != p_id && (i == 0 || ngh != candidates[i - 1].first)) {
          new_nghs[offset] = ngh;
          offset++;
        }
      }
      return;
    }
    // Sort the candidate set in reverse order according to distance from p.
    parlay::sort_inplace(candidates, less);

    size_t num_new = 0;

    size_t candidate_idx = 0;

    while (num_new < maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      node_id p_star = candidates[candidate_idx].first;
      candidate_idx++;
      if (p_star == p_id || p_star == std::numeric_limits<node_id>::max()) {
        continue;
      }

      new_nghs[num_new] = p_star;
      num_new++;

      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        node_id p_prime = candidates[i].first;
        if (p_prime != std::numeric_limits<node_id>::max()) {
          float dist_starprime =
              D->distance(data_store.find(p_star).begin(),
                          data_store.find(p_prime).begin(), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = std::numeric_limits<node_id>::max();
          }
        }
      }
    }
  }

  // robustPrune routine as found in DiskANN paper.
  // The new candidate set is added to the supplied array (new_nghs).
  void robustPrune(Graph& GG, tvec_point* p, node_id p_id,
                   parlay::sequence<pid> candidates, double alpha,
                   parlay::slice<node_id*, node_id*> new_nghs,
                   bool add = true) {
    auto less = [&](pid a, pid b) { return a.second < b.second; };
    // add out neighbors of p to the candidate set.
    if (add) {
      auto p_vertex = GG.get_vertex(p_id);
      auto map_f = [&](node_id u, node_id neighbor_id, empty_weight weight) {
        candidates.push_back(std::make_pair(
            neighbor_id, D->distance(v[neighbor_id]->coordinates.begin(),
                                     p->coordinates.begin(), d)));
        return true;
      };
      p_vertex.out_neighbors().foreach_cond(map_f);
    }
    if (candidates.size() <= maxDeg) {
      uint32_t offset = 0;
      parlay::sort_inplace(candidates, less);
      for (size_t i = 0; i < candidates.size(); ++i) {
        node_id ngh = candidates[i].first;
        if (ngh != p_id && (i == 0 || ngh != candidates[i - 1].first)) {
          new_nghs[offset] = ngh;
          offset++;
        }
      }
      return;
    }
    // Sort the candidate set in reverse order according to distance from p.
    parlay::sort_inplace(candidates, less);

    size_t num_new = 0;

    size_t candidate_idx = 0;

    while (num_new < maxDeg && candidate_idx < candidates.size()) {
      // Don't need to do modifications.
      node_id p_star = candidates[candidate_idx].first;
      candidate_idx++;
      if (p_star == p->id || p_star == std::numeric_limits<node_id>::max()) {
        continue;
      }

      new_nghs[num_new] = p_star;
      num_new++;

      for (size_t i = candidate_idx; i < candidates.size(); i++) {
        node_id p_prime = candidates[i].first;
        if (p_prime != std::numeric_limits<node_id>::max()) {
          float dist_starprime =
              D->distance(v[p_star]->coordinates.begin(),
                          v[p_prime]->coordinates.begin(), d);
          float dist_pprime = candidates[i].second;
          if (alpha * dist_starprime <= dist_pprime) {
            candidates[i].first = std::numeric_limits<node_id>::max();
          }
        }
      }
    }
  }

  // wrapper to allow calling robustPrune on a sequence of candidates
  // that do not come with precomputed distances
  void robustPrune(Graph& GG, Tvec_point<T>* p, node_id p_id,
                   parlay::sequence<node_id>& candidates, double alpha,
                   parlay::slice<node_id*, node_id*> new_nghs,
                   bool add = true) {
    parlay::sequence<pid> cc;
    auto p_vertex = GG.get_vertex(p_id);
    node_id p_ngh_size = p_vertex.out_degree();

    cc.reserve(candidates.size() + p_ngh_size);
    for (size_t i = 0; i < candidates.size(); ++i) {
      cc.push_back(std::make_pair(
          candidates[i], D->distance(v[candidates[i]]->coordinates.begin(),
                                     p->coordinates.begin(), d)));
    }
    robustPrune(GG, p, p_id, std::move(cc), alpha, new_nghs, add);
  }

    // wrapper to allow calling robustPrune on a sequence of candidates
  // that do not come with precomputed distances
  void robustPrune(Graph& GG, T* p_coords, node_id p_id,
                   parlay::sequence<node_id>& candidates, double alpha,
                   parlay::slice<node_id*, node_id*> new_nghs,
                   bool add = true) {
    parlay::sequence<pid> cc;
    auto p_vertex = GG.get_vertex(p_id);
    node_id p_ngh_size = p_vertex.out_degree();

    cc.reserve(candidates.size() + p_ngh_size);
    for (size_t i = 0; i < candidates.size(); ++i) {
      cc.push_back(std::make_pair(
          candidates[i], D->distance(data_store.find(candidates[i]).begin(),
                                     p_coords, d)));
    }
    robustPrune(GG, p_coords, p_id, std::move(cc), alpha, new_nghs, add);
  }

  // special size function
  template <class G>
  static int size_of(parlay::slice<G*, G*> nbh) {
    int size = 0;
    size_t i = 0;
    while (nbh[i] != std::numeric_limits<node_id>::max() && i < nbh.size()) {
      size++;
      i++;
    }
    return size;
  }

  void batch_insert_inplace(Graph& G, parlay::sequence<node_id>& inserts,
                            double base = 2, double max_fraction = .02,
                            bool random_order = true) {
    size_t n = v.size();
    size_t m = inserts.size();
    size_t max_batch_size = static_cast<size_t>(std::ceil(max_fraction * n));
    parlay::sequence<node_id> rperm;
    if (random_order) {
      rperm = parlay::random_permutation<node_id>(static_cast<node_id>(m),
                                                  time(NULL));
    } else {
      rperm = parlay::tabulate(m, [&](node_id i) { return i; });
    }
    auto shuffled_inserts =
        parlay::tabulate(m, [&](size_t i) { return inserts[rperm[i]]; });

    size_t floor = 0, ceiling = 0;
    uint32_t cnt_batch = 0;
    double frac=0;
    double progress_inc=.1;
    while (ceiling < m) {
      cnt_batch++;
      floor = ceiling;
      ceiling = std::min({m, size_t(floor * base) + 1, floor + max_batch_size});

      auto ind = frac*n;
      if(floor <= ind && ceiling > ind){
        frac += progress_inc;
        std::cout << "Index build " << 100*frac << "% complete" << std::endl;
      }
      // std::cout << "Start of batch insertion round: ceiling = " << ceiling
      //           << " floor = " << floor << std::endl;

      parlay::sequence<node_id> new_out = parlay::sequence<node_id>(
          maxDeg * (ceiling - floor), std::numeric_limits<node_id>::max());

      // search for each node starting from the medoid, then call
      // robustPrune with the visited list as its candidate set
      parlay::parallel_for(floor, ceiling, [&](size_t i) {
        size_t index = shuffled_inserts[i];
        T* coords = data_store.find(index).begin();
        parlay::sequence<pid> visited =
            (beam_search_2(G, coords, beamSize)).second;
        auto new_nghs =
            parlay::make_slice(new_out.begin() + maxDeg * (i - floor),
                               new_out.begin() + maxDeg * (i + 1 - floor));
        robustPrune(G, coords, index, visited, alpha, new_nghs);
      });
      // New neighbors of each point written into new_nbhs (above).
      // Not yet added to the graph G.

      // make each edge bidirectional by first adding each new edge
      // (i,j) to a sequence, then semisorting the sequence by key values
      auto to_flatten = parlay::tabulate(ceiling - floor, [&](size_t i) {
        node_id index = shuffled_inserts[i + floor];

        auto new_nghs = parlay::make_slice(new_out.begin() + maxDeg * i,
                                           new_out.begin() + maxDeg * (i + 1));

        auto edges = parlay::tabulate(size_of(new_nghs), [&](size_t j) {
          return std::make_pair((new_nghs)[j], index);
        });
        return edges;
      });

      auto KVs = parlay::tabulate(ceiling - floor, [&](size_t i) {
        node_id index = shuffled_inserts[i + floor];
        auto new_nghs = parlay::make_slice(new_out.begin() + maxDeg * i,
                                           new_out.begin() + maxDeg * (i + 1));
        size_t nghs_size = size_of(new_nghs);
        auto begin =
            (std::tuple<node_id, empty_weight>*)(new_out.begin() + maxDeg * i);
        auto tree = edge_tree(begin, begin + nghs_size);
        auto tree_ptr = tree.root;
        tree.root = nullptr;
        return std::make_tuple(index, tree_ptr);
      });

      G.insert_vertices_batch(KVs.size(), KVs.begin());
      // std::cout << "After inserts, G.num_vertices() (max node_id) = "
      //           << G.num_vertices() << std::endl;

      // TODO: update the code below:
      auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
      auto reverse_KVs =
          parlay::sequence<std::tuple<node_id, edge_node*>>(grouped_by.size());

      // std::cout << "Number of in-neighbors: " << grouped_by.size() <<
      // std::endl; finally, add the bidirectional edges; if they do not make
      // the vertex exceed the degree bound, just add them to out_nbhs;
      // otherwise, use robustPrune on the vertex with user-specified alpha
      parlay::parallel_for(0, grouped_by.size(), [&](size_t j) {
        // std::cout << "here1" << std::endl;
        auto [index, candidates] = grouped_by[j];
        // TODO: simpler case when newsize <= maxDeg.
        parlay::sequence<node_id> new_out_2(
            maxDeg, std::numeric_limits<node_id>::max());
        auto output_slice =
            parlay::make_slice(new_out_2.begin(), new_out_2.begin() + maxDeg);
        T* coords = data_store.find(index).begin();
        robustPrune(G, coords, index, candidates, alpha, output_slice);
        size_t deg = size_of(output_slice);
        auto begin = (std::tuple<node_id, empty_weight>*)new_out_2.begin();
        auto tree = edge_tree(begin, begin + deg);
        reverse_KVs[j] = {index, tree.root};
        tree.root = nullptr;
      });

      // std::cout << "ReverseKVs.size = " << reverse_KVs.size() << std::endl;
      G.insert_vertices_batch(reverse_KVs.size(), reverse_KVs.begin());
    }
  }

  Graph batch_insert_functional(Graph& G, parlay::sequence<node_id>& inserts) {
    size_t n = inserts.size();
    size_t floor = 0, ceiling = n;

    std::cout << "Functional insert size = " << n << std::endl;

    parlay::sequence<node_id> new_out = parlay::sequence<node_id>(
        maxDeg * (ceiling - floor), std::numeric_limits<node_id>::max());

    // search for each node starting from the medoid, then call
    // robustPrune with the visited list as its candidate set
    parlay::parallel_for(floor, ceiling, [&](size_t i) {
      size_t index = inserts[i];
      parlay::sequence<pid> visited =
          (beam_search(G, v[index]->coordinates.begin(), beamSize)).second;
      auto new_nghs =
          parlay::make_slice(new_out.begin() + maxDeg * (i - floor),
                             new_out.begin() + maxDeg * (i + 1 - floor));
      robustPrune(G, v[index], index, visited, alpha, new_nghs);
    });
    // New neighbors of each point written into new_nbhs (above).
    // Not yet added to the graph G.

    // make each edge bidirectional by first adding each new edge
    // (i,j) to a sequence, then semisorting the sequence by key values
    auto to_flatten = parlay::tabulate(ceiling - floor, [&](size_t i) {
      node_id index = inserts[i + floor];

      auto new_nghs = parlay::make_slice(new_out.begin() + maxDeg * i,
                                         new_out.begin() + maxDeg * (i + 1));

      auto edges = parlay::tabulate(size_of(new_nghs), [&](size_t j) {
        return std::make_pair((new_nghs)[j], index);
      });
      return edges;
    });

    auto KVs = parlay::tabulate(ceiling - floor, [&](size_t i) {
      node_id index = inserts[i + floor];
      auto new_nghs = parlay::make_slice(new_out.begin() + maxDeg * i,
                                         new_out.begin() + maxDeg * (i + 1));
      size_t nghs_size = size_of(new_nghs);
      auto begin =
          (std::tuple<node_id, empty_weight>*)(new_out.begin() + maxDeg * i);
      auto tree = edge_tree(begin, begin + nghs_size);
      auto tree_ptr = tree.root;
      tree.root = nullptr;

      return std::make_tuple(index, tree_ptr);
    });

    Graph new_G = G.insert_vertices_batch_functional(KVs.size(), KVs.begin());
    // std::cout << "After inserts, G.num_vertices() (max node_id) = "
    //           << new_G.num_vertices() << std::endl;

    // TODO: update the code below:
    auto grouped_by = parlay::group_by_key(parlay::flatten(to_flatten));
    auto reverse_KVs =
        parlay::sequence<std::tuple<node_id, edge_node*>>(grouped_by.size());

    // std::cout << "Number of in-neighbors: " << grouped_by.size() <<
    // std::endl; finally, add the bidirectional edges; if they do not make the
    // vertex exceed the degree bound, just add them to out_nbhs; otherwise, use
    // robustPrune on the vertex with user-specified alpha
    parlay::parallel_for(0, grouped_by.size(), [&](size_t j) {
      auto [index, candidates] = grouped_by[j];
      // TODO: simpler case when newsize <= maxDeg.
      parlay::sequence<node_id> new_out_2(maxDeg,
                                          std::numeric_limits<node_id>::max());
      auto output_slice =
          parlay::make_slice(new_out_2.begin(), new_out_2.begin() + maxDeg);
      robustPrune(new_G, v[index], index, candidates, alpha, output_slice);
      size_t deg = size_of(output_slice);
      auto begin = (std::tuple<node_id, empty_weight>*)new_out_2.begin();
      auto tree = edge_tree(begin, begin + deg);
      reverse_KVs[j] = {index, tree.root};
      tree.root = nullptr;
    });

    std::cout << "ReverseKVs.size = " << reverse_KVs.size() << std::endl;
    //    Graph newer_G =
    //    new_G.insert_vertices_batch_functional(reverse_KVs.size(),
    //                                                           reverse_KVs.begin());
    //    return newer_G;

    new_G.insert_vertices_batch(reverse_KVs.size(), reverse_KVs.begin());
    return new_G;
  }

  parlay::sequence<float> centroid_helper(slice_tvec a) {
    if (a.size() == 1) {
      parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
      for (unsigned i = 0; i < d; i++)
        centroid_coords[i] = static_cast<float>((a[0]->coordinates)[i]);
      return centroid_coords;
    } else {
      size_t n = a.size();
      parlay::sequence<float> c1;
      parlay::sequence<float> c2;
      parlay::par_do_if(
          n > 1000, [&]() { c1 = centroid_helper(a.cut(0, n / 2)); },
          [&]() { c2 = centroid_helper(a.cut(n / 2, n)); });
      parlay::sequence<float> centroid_coords = parlay::sequence<float>(d);
      for (unsigned i = 0; i < d; i++) {
        float result = (c1[i] + c2[i]) / 2;
        centroid_coords[i] = result;
      }
      return centroid_coords;
    }
  }

  // TODO(laxmand): Why fvec point here?
  tvec_point* medoid_helper(tvec_point* centroid, slice_tvec a) {
    if (a.size() == 1) {
      return a[0];
    } else {
      size_t n = a.size();
      tvec_point* a1;
      tvec_point* a2;
      parlay::par_do_if(
          n > 1000, [&]() { a1 = medoid_helper(centroid, a.cut(0, n / 2)); },
          [&]() { a2 = medoid_helper(centroid, a.cut(n / 2, n)); });
      float d1 = D->distance(centroid->coordinates.begin(),
                             a1->coordinates.begin(), d);
      float d2 = D->distance(centroid->coordinates.begin(),
                             a2->coordinates.begin(), d);
      if (d1 < d2)
        return a1;
      else
        return a2;
    }
  }

  // computes the centroid and then assigns the approx medoid as the point in v
  // closest to the centroid
  node_id find_approx_medoid() {
    parlay::sequence<float> centroid = centroid_helper(parlay::make_slice(v));
    auto c = parlay::tabulate(
        centroid.size(), [&](size_t i) { return static_cast<T>(centroid[i]); });
    tvec_point centroidp = tvec_point();
    centroidp.coordinates = parlay::make_slice(c);
    medoid = medoid_helper(&centroidp, parlay::make_slice(v));
    std::cout << "Medoid ID: " << medoid->id << std::endl;
    return medoid->id;
  }

};
