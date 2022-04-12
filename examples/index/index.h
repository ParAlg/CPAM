#pragma once

#include <vector>

#include <cpam/cpam.h>
#include <pam/pam.h>
// using namespace cpam;
#include "index_encoder.h"
#include <pam/get_time.h>
#include <pam/parse_command_line.h>

//#include <pam/pam.h>
//#include <pam/get_time.h>
//#include <pam/parse_command_line.h>

auto add = [] (float a, float b) -> float {
    return a + b;};

using token = parlay::sequence<char>;

struct inv_index {
  using doc_id = int;
  using weight = int;

  using post_elt = pair<doc_id, weight>;

  struct doc_entry {
    using key_t = doc_id;
    using val_t = weight;
    static inline bool comp(key_t a, key_t b) { return a < b;}
    using aug_t = weight;
    static aug_t get_empty() {return 0;}
    static aug_t from_entry(key_t k, val_t v) {return v;}
    static aug_t combine(aug_t a, aug_t b) {
      return (b > a) ? b : a;}
    using entry_t = std::pair<key_t,val_t>;
  };

  //using post_list = aug_map<doc_entry>;
  //using post_list = aug_map<doc_entry>;
#ifdef USE_DIFF_ENCODING
  using post_list = cpam::aug_map<doc_entry, 16, cpam::diffencoded_index_encoder>;
#else
#ifdef USE_PAM
  using post_list = aug_map<doc_entry>;
#else
  using post_list = cpam::aug_map<doc_entry, 16>;
#endif
#endif
  using post_list_node = typename post_list::node;

  struct token_entry {
    using key_t = token;
    using val_t = post_list;
    static inline bool comp(const key_t& a, const key_t& b) {
      size_t m = std::min(a.size(), b.size());
      for (size_t i = 0; i < m; i++)
	if (a[i] != b[i]) return a[i] < b[i];
      return a.size() < b.size();
    }
  };

  using index_elt = pair<token, post_elt>;
  using index_list = pair<token, post_list>;

#ifdef USE_PAM
  using index = pam_map<token_entry>;
#else
  using index = cpam::pam_map<token_entry, 2>;
#endif
  using index_node = typename index::node;

  index idx;

  inv_index(parlay::sequence<index_elt> const &S) {
    timer t("build index", false);
    //size_t n = S.size();
    auto reduce = [&] (parlay::slice<post_elt*,post_elt*> R) {
      // should be optimized to build sequentially in place
      //return post_list(R, add, true); };
      return post_list(R, add); };
    idx = index::multi_insert_reduce(index(), S, reduce);
    t.next("build");
  }

  inv_index() {}

  post_list get_list(const token w) {
    std::optional<post_list> p = idx.find(w);
    if (p) return *p;
    else return post_list();
  }

  static post_list And(post_list a, post_list b) {
    return post_list::map_intersect(std::move(a),std::move(b),add);}

  static post_list Or(post_list a, post_list b) {
    return post_list::map_union(a,b,add);}

  static post_list And_Not(post_list a, post_list b) {
    return post_list::map_difference(a,b);}

  void print_index_size() {
#ifdef USE_PAM
    size_t num_outer = 0;
    size_t num_inner = 0;
    auto fn = [&] (const auto& et) {
      num_outer++;
      num_inner += std::get<1>(et).size();
    };
    idx.foreach_seq(idx, fn);
    size_t total_size = num_outer*sizeof(index_node) + num_inner*sizeof(post_list_node);
    std::cout << "Num outer nodes = " << num_outer << " Num inner nodes = " << num_inner << std::endl;
    std::cout << "Total size in bytes: " << total_size << std::endl;
#else
    auto fn_noop = [&] (const auto& et) {
      return 0;
    };
    size_t outer_bytes = idx.size_in_bytes(fn_noop);
    std::cout << "Outer bytes = " << outer_bytes << std::endl;

    auto fn_inner = [&] (const auto& et) {
      auto inner_noop = [&] (const auto& et) {
        return 0;
      };
      size_t inner_size = std::get<1>(et).size_in_bytes(inner_noop);
      return inner_size;
    };
    size_t inner_bytes = idx.size_in_bytes(fn_inner);
    std::cout << "Inner bytes = " << inner_bytes << std::endl;
    std::cout << "Total bytes = " << (outer_bytes + inner_bytes) << std::endl;
#endif
  }

  size_t index_size() {
#ifdef USE_PAM
    size_t num_outer = 0;
    size_t num_inner = 0;
    auto fn = [&] (const auto& et) {
      num_outer++;
      num_inner += std::get<1>(et).size();
    };
    idx.foreach_seq(idx, fn);
    size_t total_size = num_outer*sizeof(index_node) + num_inner*sizeof(post_list_node);
    return total_size;
#else
    auto fn_noop = [&] (const auto& et) {
      return 0;
    };
    size_t outer_bytes = idx.size_in_bytes(fn_noop);

    auto fn_inner = [&] (const auto& et) {
      auto inner_noop = [&] (const auto& et) {
        return 0;
      };
      size_t inner_size = std::get<1>(et).size_in_bytes(inner_noop);
      return inner_size;
    };
    size_t inner_bytes = idx.size_in_bytes(fn_inner);
    return outer_bytes + inner_bytes;
#endif
  }

  vector<post_elt> top_k(const post_list& a, int k) {
    int l = std::min<int>(k,a.size());
    vector<post_elt> vec(l);
    post_list b = a;
//    std::cout << "Starting top_k, size = " << a.size() << std::endl;
//    std::cout << "ref_cnt is now: " << b.ref_cnt() << std::endl;
    for (int i=0; i < l; i++) {
      weight m = b.aug_val();
      auto f = [m] (weight v) {return v < m;};
      auto [u, v] = *b.aug_select(f);
      vec[i] = post_elt(u, v);
//      std::cout << "key = " << vec[i].first << " b.size = " << b.size() << "ref_cnt = " << b.ref_cnt() << std::endl;
      b = post_list::remove(std::move(b),vec[i].first);
 //     std::cout << "After remove, b.size is now: " << b.size() << " ref_cnt = " << b.ref_cnt() << std::endl;
 //     std::cout << "is_compressed:" << b.root_is_compressed() << std::endl;
    }
    return vec;
  }

};
