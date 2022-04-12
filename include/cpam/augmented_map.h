#pragma once

using namespace std;

// *******************************************
//   AUG MAPS
// *******************************************

namespace cpam {

template <class _Entry, class Join_Tree>
struct aug_map_ : private map_<_Entry, Join_Tree> {
public:
  using Map = map_<_Entry, Join_Tree>;
  using Entry = typename Map::Entry;
  using Tree = augmented_ops<typename Map::Tree>;
  using node = typename Tree::node;
  using E = typename Entry::entry_t;
  using K = typename Entry::key_t;
  using V = typename Entry::val_t;
  using A = typename Entry::aug_t;
  using M = aug_map_;
  using GC = typename Map::GC;
  using maybe_V = std::optional<V>;
  using maybe_E = std::optional<E>;
  using ptr = typename GC::ptr;
  using Map::kNodeLimit;

  template<class F>
  static M aug_filter(M m, const F& f) {
    return M(Tree::aug_filter(m.get_root(), f)); }

  // extract the augmented values
  A aug_val() const { return Tree::aug_val(Map::root); }

  A aug_left (const K& key) {
    typename Tree::aug_sum_t a;
    Tree::aug_sum_left(Map::root, key, a);  // NOT using ptr.
    return a.result;}

  A aug_right(const K& key) {
    typename Tree::aug_sum_t a;
    Tree::aug_sum_right(Map::root, key, a);
    return a.result;}

  A aug_range(const K& key_left, const K& key_right) {
    typename Tree::aug_sum_t a;
    Tree::aug_sum_range(Map::root, key_left, key_right, a);
    return a.result;}

  // just side effecting
  template <class restricted_sum>
  void range_sum(const K& key_left,
		 const K& key_right,
		 restricted_sum& rs) {
    Tree::aug_sum_range(Map::root, key_left, key_right, rs);
  }

  template <class Func>
  maybe_E aug_select(Func f) {
    return Tree::aug_select(Map::root, f); }

  static M insert_lazy(M m, const E& p) {
    auto replace = [] (const V& a, const V& b) {return b;};
    return M(Tree::insert_lazy(m.get_root(), p, replace)); }

  // for coercing a map to an aug_map, should be a better way
  static M to_aug(Map&& m) {
    aug_map_ x;
    x.root = m.root; m.root = NULL; return x;}
  aug_map_() : Map() { }
  // install Map's constructors
  using Map::Map;

  template <class Func>
  static M insert(M m, const E& p, const Func& f) {
    return to_aug(Map::insert(std::move(m),p,f));}
  static M insert(M m, const E& p) {//cout << "ins a " << endl;
    return to_aug(Map::insert(std::move(m),p));}
  static M remove(M m, const K& k) {return to_aug(Map::remove(std::move(m), k));}
  template<class F>
  static M filter(M m, const F& f) {return to_aug(Map::filter(std::move(m), f));}
  static M multi_insert(M m, parlay::sequence<E> const &SS) {
    return to_aug(Map::multi_insert(std::move(m), SS));}
  template<class Seq, class BinOp>
  static M multi_insert_sorted(M m, Seq const &SS, BinOp f) {
    return to_aug(Map::multi_insert_sorted(std::move(m), SS, f));
  }
  template<class Seq, class CombineOp, class MapOp>
  static M multi_insert_sorted_map(M m, Seq &SS, CombineOp combine_op, MapOp map_op) {
    return to_aug(Map::multi_insert_sorted_map(std::move(m), SS, combine_op, map_op));
  }
  template<class Seq, class CombineOp>
  static M multi_delete_sorted_map(M m, Seq &SS, CombineOp combine_op) {
    return to_aug(Map::multi_delete_sorted_map(std::move(m), SS, combine_op));
  }
  // delete multiple keys from an array
  template<class Seq>
  static M multi_delete_sorted(M m, Seq const &SS) {
    return to_aug(Map::multi_delete_sorted(std::move(m), SS));
  }

  // update multiple entries from a sorted array
  template<class Seq, class Bin_Op>
  static bool multi_update_sorted_inplace(M& m, Seq &SS, Bin_Op f) {
    return Tree::multi_update_sorted_inplace(m.root, SS.begin(), SS.size(), f);
  }

  template <class Seq>
  static parlay::sequence<V> multi_find_sorted(M& m, Seq const &SS) {
    return Map::multi_find_sorted(m, SS);
  }

  template<class Bin_Op>
  static M multi_insert_combine(M m, parlay::sequence<E> S, Bin_Op f,  // ?? should it be &
				bool seq_inplace = false) {
    return to_aug(Map::multi_insert_combine(std::move(m), S, f, seq_inplace));}
  template<class Val, class Reduce>
  static M multi_insert_reduce(M m, const parlay::sequence<pair<K,Val>>& S, Reduce g) {  // ?? should it be &
    return to_aug(Map::multi_insert_reduce(std::move(m), S, g)); }

  template<class M1, class M2, class F>
  static M map_intersect(M1 a, M2 b, const F& op) {
    return to_aug(Map::map_intersect(std::move(a), std::move(b), op));}
  static M map_intersect(M a, M b) {return to_aug(Map::map_intersect(std::move(a), std::move(b)));}
  template<class F>
  static M map_union(M a, M b, const F& op) {return to_aug(Map::map_union(std::move(a), std::move(b), op));}
  static M map_union(M a, M b) {return to_aug(Map::map_union(std::move(a), std::move(b)));}
  static M map_difference(M a, M b) {return to_aug(Map::map_difference(std::move(a), std::move(b)));}

  static M range(M& a, K kl, K kr) {return to_aug(Map::range(a,kl,kr));}
  template<class Ma, class F>
  static M map(Ma& a, const F f) {return to_aug(Map::map(a, f));}
  static M subseq(M& a, size_t left, size_t right) {return to_aug(Map::subseq(a,left,right));}

  static size_t size(node* r) { return Map::size(r); }

  static void entries(M m, E* out) { Map::entries(std::move(m),out);}
  static parlay::sequence<E> entries(M m, size_t granularity=kNodeLimit) { return Map::entries(std::move(m));}
  template <class outItter>
  static void keys(M m, outItter out) {Map::keys(std::move(m),out);}
  static void keys_to_array(M m, K* out) {Map::keys_to_array(std::move(m),out);}
  static parlay::sequence<K> keys(M m, size_t granularity=kNodeLimit) {
	  return Map::keys(m, granularity);
  }
  bool operator == (const M& m) { return Map::operator==(m);}

  template<class R, class F>
  static typename R::T map_reduce(const M& m, const F& f, const R& r,
				  size_t grain=kNodeLimit) {
    return Map::template map_reduce<R>(m, f, r, grain);}
  template<class F>
  static void map_index(M m, const F& f, size_t granularity = kNodeLimit,
			size_t start=0) {
    Map::map_index(m, f, granularity, start);
  }
  template<class Ma, class F>
  static M map_set(Ma a, const F& f) {return to_aug(Map::map_set(a, f));}

  template <class F>
  static void foreach_index(const M& m, const F& f, size_t start=0,
			    size_t granularity = kNodeLimit) {
    Map::foreach_index(m, f, start, granularity); }
  template <class F>
  static void foreach_index_2(const M& m, const F& f) {
    Tree::foreach_index_2(m.root, f);
  }
  template <class F>
  static void foreach_cond(M m, F f, size_t start=0,
			    size_t granularity = kNodeLimit) {
    Map::foreach_cond(m, f, start, granularity); }
  template <class F, class C>
  static void foreach_cond_par(M& m, F f, C c) {
    Map::foreach_cond_par(m, f, c); }
  // apply function f on all entries sequentially
  template <class F>
  static void foreach_seq(const M& m, const F& f) {
    Map::foreach_seq(m, f);
  }
public:
  using Map::from_sorted;
  using Map::size;
  using Map::is_empty;
  using Map::init;
  using Map::reserve;
  using Map::finish;
  using Map::clear;
  using Map::find;
  using Map::contains;
  using Map::rank;
  using Map::select;
  using Map::subseq;
  using Map::root;
  using Map::get_root;
  using Map::update;
  using Map::insert;
  using Map::check_balance;
  using Map::root_is_compressed;
  using Map::check_structure;
  using Map::size_in_bytes;
  using Map::node_stats;
  using Map::iterate_seq;
  using Map::ref_cnt;
};

// creates a key-value pair for the entry, and redefines from_entry
template <class entry>
struct aug_map_full_entry : entry {
  using val_t = typename entry::val_t;
  using key_t = typename entry::key_t;
  using aug_t = typename entry::aug_t;
  using entry_t = std::tuple<key_t,val_t>;
  static inline key_t get_key(const entry_t& e) {return std::get<0>(e);}
  static inline val_t get_val(const entry_t& e) {return std::get<1>(e);}
  static inline void set_val(entry_t& e, const val_t& v) {std::get<1>(e) = v;}
  static inline entry_t to_entry(const key_t& k, const val_t& v) {
    return std::make_tuple(k, v);};
  static inline aug_t from_entry(const entry_t& e) {
    return entry::from_entry(std::get<0>(e), std::get<1>(e));}
};

template <class _Entry, size_t BlockSize=128, class Encoder=default_entry_encoder, class Balance = weight_balanced_tree>
using aug_map =
    aug_map_<aug_map_full_entry<_Entry>,
             typename Balance::template balance<
                 aug_node<typename Balance::data,
                          aug_map_full_entry<_Entry>,
                          typename Encoder::template encoder<aug_map_full_entry<_Entry>, /* is_aug = */ true>,
                          BlockSize
                        >>>;

template <class _Entry, size_t BlockSize=256, class Balance=weight_balanced_tree>
using diff_encoded_aug_map = aug_map<_Entry, BlockSize, diffencoded_entry_encoder, Balance>;

// creates a key-value pair for the entry, and redefines from_entry
template <class entry>
struct aug_set_full_entry : entry {
  using val_t = bool;  // not used
  using key_t = typename entry::key_t;
  using aug_t = typename entry::aug_t;
  using entry_t = key_t;
  static inline key_t get_key(const entry_t& e) { return e; }
  static inline val_t get_val(const entry_t& e) { return 0; }
  static inline void set_val(entry_t& e, const val_t& v) {}
};

// _Entry needs:
//    key_t, aug_t,
//    comp(key_t, key_t) -> bool,
//    from_entry(key_t) -> aug_t,
//    get_empty() -> aug_tm,
//    combine(aug_t, aug_t) -> aug_t
template <class _Entry, size_t BlockSize=256, class Encoder=default_entry_encoder, class Balance = weight_balanced_tree>
using aug_set =
    aug_map_<aug_set_full_entry<_Entry>,
             typename Balance::template balance<aug_node<
                 typename Balance::data,
                 aug_set_full_entry<_Entry>,
                 typename Encoder::template encoder<aug_set_full_entry<_Entry>, /* is_aug = */ true>,
                 BlockSize
             >>>;

}  // namespace cpam
