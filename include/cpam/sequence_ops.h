#pragma once
#include "gc.h"
#include "utils.h"

// *******************************************
//   SEQUENCES
// *******************************************

namespace cpam {

template<class Tree>
struct sequence_ops : Tree {
  using node = typename Tree::node;
  using regular_node = typename Tree::regular_node;
  using compressed_node = typename Tree::compressed_node;
  using ET = typename Tree::ET;
  using GC = gc<Tree>;
  using ptr = typename GC::ptr;
  using Tree::B;
  using Tree::kBaseCaseSize;
  using Tree::kNodeLimit;

  using expose_tuple = std::tuple<ptr, ET&, ptr, regular_node*>;
  using expose_simple_tuple = std::tuple<node*, ET&, node*>;

  // Takes a compressed node and returns the node's elements as a tree
  static regular_node* to_tree_impl(ET* A, uint32_t n, uint32_t depth = 0) {
    if (n <= 0) return (regular_node*)Tree::empty();
    if (n == 1) return Tree::single(A[0]);
    if (n == B && depth > 0) {
      return (regular_node*)(Tree::make_single_compressed_node(A, n));
    }

    size_t mid = n/2;
    regular_node* m = Tree::make_regular_node(A[mid]);

    size_t l_size = mid, r_size = n-mid-1;
    auto P = utils::fork<regular_node*>(n >= kNodeLimit,
      [&]() {return to_tree_impl(A, l_size, depth + 1);},
      [&]() {return to_tree_impl(A+mid+1, r_size, depth + 1);});
    return Tree::regular_node_join(P.first, P.second, m);
  }

  static regular_node* to_tree(compressed_node* c) {
    ET tmp_arr[2*B];
    ET* arr = Tree::compressed_node_elms(c, tmp_arr);
    size_t arr_size = c->s;
    assert(arr_size > 0);
    return to_tree_impl(arr, arr_size);
  }

  static expose_simple_tuple expose_simple(node* p) {
    if (Tree::is_regular(p)) {
      regular_node* rp = (regular_node*)p;
      return expose_simple_tuple(rp->lc, Tree::get_entry(p), rp->rc);
    } else {
      std::cout << "error in expose_simple: applied to a compressed node" << std::endl;
      return expose_simple_tuple(nullptr, Tree::get_entry(p), nullptr);
    }
  }

  static expose_tuple expose(ptr a) {
    node* p = ptr::strip_flag(a.p);
    bool extra = a.extra();
    bool other_real_refs = (Tree::ref_cnt(p) > 1);
    bool copy = extra || other_real_refs;
    bool is_regular = Tree::is_regular(p);
    if (is_regular) {
      regular_node* ret = nullptr;
      regular_node* rp = (regular_node*)p;
      if (!copy) {
        a.p = nullptr;
        ret = Tree::cast_to_regular(p);
      }
      auto ret_tup = expose_tuple(ptr(rp->lc, copy), Tree::get_entry(rp), ptr(rp->rc, copy), ret);
      return ret_tup;
    } else {
      // convert to tree
      //if (print) {
      //  std::cout << "Hit a compressed node, converting to tree. Size = " << Tree::size(p) << std::endl;
      //  size_t sz = 0;
      //  auto fn = [&] (const auto& et) {
      //    std::cout << std::get<0>(et) << " " << sz << std::endl;
      //    sz++;
      //  };
      //  Tree::iterate_seq(p, fn);
      //}
      regular_node* root = to_tree((compressed_node*)p);
      //if (print) {
      //  auto fn = [&] (const auto& et) {
      //    std::cout << et.first << std::endl;
      //  };
      //  std::cout << "printing new root" << std::endl;
      //  Tree::iterate_seq(root, fn);
      //}
      return expose_tuple(ptr(root->lc, false), Tree::get_entry(root), ptr(root->rc, false), root);
    }
  }

  static node* join(node* l, ET e, node* r, regular_node* root) {
    if (root == nullptr) {
      root = Tree::make_regular_node(e);
    }
    return Tree::node_join(l, r, root);
  }

  static node* regular_join(node* l, ET e, node* r, regular_node* root) {
    if (root == nullptr) {
      root = Tree::make_regular_node(e);
    }
    return Tree::regular_node_join(Tree::cast_to_regular(l),
        Tree::cast_to_regular(r), root);
  }


  static node_size_t depth(node* a) {
    if (a == NULL) return 0;
    auto P = utils::fork<node_size_t>(Tree::size(a) >= kNodeLimit,
		[&]() {return depth(a->lc);},
		[&]() {return depth(a->rc);});
    return std::max(P.first, P.second) + 1;
  }

  static bool check_balance(node* a) {
    if (a == NULL) return true;
    auto P = utils::fork<bool>(Tree::size(a) >= kNodeLimit,
		[&]() {return check_balance(a->lc);},
		[&]() {return check_balance(a->rc);});
    return Tree::is_balanced(a) && P.first && P.second;
  }

  // TODO: add base case.
  static std::optional<ET> select_compressed(node* a, size_t rank) {
    size_t i=0;
    std::optional<ET> ret;
    auto fn = [&] (const ET& a) -> bool {
      if (i == rank) {
        ret = a;
        return false;
      }
      i++;
      return true;
    };
    Tree::iterate_cond(a, fn);
    return ret;
  }

  // TODO: add base case.
  static std::optional<ET> select(node* a, size_t rank) {
    if (a == nullptr) return {};
    //size_t size = b.size();
    if (Tree::is_compressed(a)) return select_compressed(a, rank);
    auto an = (regular_node*)a;
    //std::cout << "Size = " << size << " e = " << e << std::endl;
    size_t left_size = Tree::size(an->lc);;
    std::optional<ET> ret;
    if (rank > left_size) {
      ret = select(an->rc, rank - left_size - 1);
    } else if (rank < left_size) {
      ret = select(an->lc, rank);
    } else {
      ret = Tree::get_entry(an);
    }
    return ret;
  }

  static node* take_compressed(node* a, size_t rank) {
    ET arr[2*B];
    size_t offset = 0;
    auto copy_f = [&] (ET e) {
      if (offset < rank) {
        parlay::move_uninitialized(arr[offset++], e);
      }
    };
    Tree::iterate_seq(a, copy_f);
    if (offset >= B) {
      return Tree::make_compressed(arr, offset);
    } else {
      return to_tree_impl(arr, offset);
    }
  }

  static node* take(node* a, size_t rank) {
    if (a == nullptr) return NULL;
    if (Tree::is_compressed(a)) return take_compressed(a, rank);

    auto ar = (regular_node*)a;
    size_t left_size = Tree::size(ar->lc);
    if (rank < left_size) {
      return take(ar->lc, rank);
    }
    regular_node* join = Tree::make_regular_node(Tree::get_entry(ar));
    node* r = take(ar->rc, rank - left_size - 1);
    GC::increment(ar->lc);
    return Tree::node_join(ar->lc, r, join);
  }

//  static node* suffix_compressed(node* a, size_t rank) {
//    ET arr[2*B];
//    size_t offset = 0;
//    size_t idx = 0;
//    auto copy_f = [&] (ET e) {
//      if (idx >= rank) {
//        parlay::move_uninitialized(arr[offset++], e);
//      }
//      idx++;
//    };
//    Tree::iterate_seq(a, copy_f);
//    if (offset >= B) {
//      return Tree::make_compressed(arr, offset);
//    } else {
//      return to_tree_impl(arr, offset);
//    }
//  }
//
//  // Suffix starting from rank
//  static node* suffix(node* a, size_t rank) {
//    if (a == nullptr) return NULL;
//    if (Tree::is_compressed(a)) return suffix_compressed(a, rank);
//
//    auto ar = (regular_node*)a;
//    size_t left_size = Tree::size(ar->lc);
//    if (rank > left_size) return suffix(ar->rc, rank - left_size - 1);
//    regular_node* join = Tree::make_regular_node(Tree::get_entry(ar));
//    node* l = suffix(ar->lc, rank);
//    GC::increment(ar->rc);
//    return Tree::node_join(l, ar->rc, join);
//  }
//
//  static node* subseq_compressed(node* a, size_t left, size_t right) {
//    ET arr[2*B];
//    size_t offset = 0;
//    size_t idx = 0;
//    auto copy_f = [&] (ET e) {
//      if ((left <= idx) && (idx < right)) {
//        parlay::move_uninitialized(arr[offset++], e);
//      }
//      idx++;
//    };
//    Tree::iterate_seq(a, copy_f);
//    if (offset >= B) {
//      return Tree::make_compressed(arr, offset);
//    } else {
//      return to_tree_impl(arr, offset);
//    }
//  }

  static std::tuple<node*, size_t, size_t> subseq_root(node* b, size_t left, size_t right) {
    if (!b) return {NULL, left, right};
    if (Tree::is_compressed(b)) {
      return {b, left, right};
    }
    auto rb = Tree::cast_to_regular(b);
    size_t left_size = Tree::size(rb->lc);
    if (right < left_size) {
      return subseq_root(rb->lc, left, right);
    } else if (left > left_size + 1) {
      return subseq_root(rb->rc, left - left_size - 1, right - left_size - 1);
    }
    return {b, left, right};
  }

//  static node* subseq(node* b, const K& low, const K& high) {
//    node* r;
//    std::tie(r, low, high) = subseq_root(b, low, high);
//    if (!r) return NULL;
//    if (Seq::is_compressed(r)) {
//      subseq_compressed(r, low, high);
//    } else {
//      auto [lc, e, rc] = Seq::expose_simple(r);
//      return Seq::join(suffix(lc, low), e, take(rc, high), nullptr);
//    }
//  }


  template <class F>
  static node* entry_bc(node* r, F& f) {
    ET stack[2*B];
    size_t offset = 0;
    auto copy_f = [&] (ET a) {
      if (f(a)) {
        parlay::move_uninitialized(stack[offset++], a);
      }
    };
    Tree::iterate_seq(r, copy_f);
    assert(offset <= 2*B);

    if (offset < B) {
      return to_tree_impl((ET*)stack, offset);
    } else {
      return Tree::make_compressed(stack, offset);
    }
  }

  // Keep elements from rank (inclusive) to end. take o rev
  static node* suffix(node* a, const size_t rank) {
    if (!a) return NULL;
    if (Tree::is_compressed(a)) {
      size_t i = 0;
      auto comp = [&] (const ET& e) {
        bool ret = i >= rank;
        i++;
        return ret; };
      return entry_bc(a, comp);
    }
    auto ar = (regular_node*)a;
    size_t left_size = Tree::size(ar->lc);
    // TODO(laxman): check edge case here (+1 or not?)
    if (rank > left_size + 1) return suffix(ar->rc, rank - left_size - 1);
    regular_node* join = Tree::make_regular_node(Tree::get_entry(ar));
    node* l = suffix(ar->lc, rank);
    GC::increment(ar->rc);
    return Tree::node_join(l, ar->rc, join);
  }

  static node* subseq(node* b, size_t low, size_t high) {
    node* r;
    std::tie(r, low, high) = subseq_root(b, low, high);
    if (!r) return NULL;
    if (Tree::is_compressed(r)) {
      // Build tree only on elements in the range in the leaf.
      size_t i=0;
      auto comp = [&] (const ET& k) {
        bool ret = (low <= i) && (i <= high);
        i++;
        return ret;
      };
      return entry_bc(r, comp);
    } else {
      auto [lc, e, rc] = expose_simple(r);
      return join(suffix(lc, low), e, take(rc, high - Tree::size(lc) - 1), nullptr);
    }
  }

  static node* join2_i(ptr b1, ptr b2) {
    if (b1.empty()) return b2.node_ptr();
    if (b2.empty()) return b1.node_ptr();

    if (b1.size() > b2.size()) {
      if (b1.is_compressed() || Tree::will_be_compressed(b1.size() + b2.size())) {
        auto ret = Tree::make_compressed(b1.node_ptr(), b2.node_ptr());
//        Tree::check_structure(ret);
        return ret;
      }
      auto[lc, e, rc, root] = expose(std::move(b1));
      node* l = lc.node_ptr();
      node* r = join2_i(std::move(rc), std::move(b2));
      auto ret = join(l, e, r, root);
//      Tree::check_structure(ret);
      return ret;
    } else {
      if (b2.is_compressed() || Tree::will_be_compressed(b1.size() + b2.size())) {
        return Tree::make_compressed(b1.node_ptr(), b2.node_ptr());
      }
      auto[lc, e, rc, root] = expose(std::move(b2));
      node* l = join2_i(std::move(b1), std::move(lc));
      assert(l != root);
      node* r = rc.node_ptr();
      auto ret = join(l, e, r, root);
//      Tree::check_structure(ret);
      return ret;
    }
  }

  static node* join2(node* b1, node* b2) {
    return join2_i(ptr(b1), ptr(b2));
  }

  template<class InTree, class Func>
  static node* map_compressed(typename InTree::ptr b, const Func& f) {
    assert(b.is_compressed());
    assert(b.size() <= 2*B);
    ET arr[2*B];
    size_t offset = 0;
    auto copy_f = [&] (typename InTree::ET e) {
      auto y = f(e);
      parlay::move_uninitialized(arr[offset++], y);
    };
    Tree::iterate_seq(b.node_ptr(), copy_f);
    return Tree::make_compressed(arr, offset);
  }

  template<class InTree, class Func>
  static node* map(typename InTree::ptr b, const Func& f) {
    if (b.empty()) return NULL;
    if (b.is_compressed()) {
      return map_compressed<InTree, Func>(std::move(b), f);
    }
    // size_t n = b.size();
    auto [lc, e, rc, m] = expose(std::move(b));
    auto P = utils::fork<node*>(true,
       [&] () {return map<InTree>(std::move(lc), f);},
       [&] () {return map<InTree>(std::move(rc), f);});
    auto y = f(e);
    regular_node* r = Tree::single(y);
    return Tree::node_join(P.first, P.second, r);
  }


  // F applied to entries.
  template <class F>
  static size_t size_in_bytes(node* a, const F& f) {
    if (!a) return 0;

    size_t total = 0;
    if (Tree::is_compressed(a)) {
      auto c = Tree::cast_to_compressed(a);
      total += c->size_in_bytes;
      auto fn = [&] (const auto& et) {
        total += f(et);
      };
      Tree::iterate_seq(c, fn);
    } else {
      size_t n = Tree::size(a);
      auto r = Tree::cast_to_regular(a);
      auto P = utils::fork<size_t>(n >= kNodeLimit,
        [&] () {return size_in_bytes(r->lc, f);},
        [&] () {return size_in_bytes(r->rc, f);});

      total += P.first + P.second;
      total += sizeof(regular_node);
      total += f(Tree::get_entry(r));
    }
    return total;
  }

  // Returns a pair of the (number of internal nodes, number of leaf nodes,
  // total size of leaf nodes).
  static std::tuple<size_t, size_t, size_t> node_stats(node* a) {
    if (!a) return {0, 0, 0};

    if (Tree::is_compressed(a)) {
      return {0, 1, Tree::size(a)};
    }
    size_t internal = 0; size_t leaf = 0; size_t leaf_sizes = 0;
    size_t n = Tree::size(a);
    auto r = Tree::cast_to_regular(a);
    auto P = utils::fork<std::tuple<size_t, size_t, size_t>>(n >= kNodeLimit,
      [&] () {return node_stats(r->lc);},
      [&] () {return node_stats(r->rc);});

    internal += std::get<0>(P.first) + std::get<0>(P.second) + 1;
    leaf += std::get<1>(P.first) + std::get<1>(P.second);
    leaf_sizes += std::get<2>(P.first) + std::get<2>(P.second);

    return {internal, leaf, leaf_sizes};
  }

  template <class F>
  static std::tuple<bool, ET, ET> is_sorted_impl(node* a, const F& less) {
    if (a == nullptr) return {true, ET(), ET()};
    if (Tree::is_compressed(a)) {
      ET left, right;
      size_t i = 0;
      bool sorted = true;
      auto fn = [&] (const ET& a) {
        if (i == 0) { left = a; }
        if (i > 0) { sorted &= less(right, a); }
        right = a;
        i++;
      };
      Tree::iterate_seq(a, fn);
      return {sorted, left, right};
    }

    auto an = (regular_node*)a;
    auto P = utils::fork<std::tuple<bool, ET, ET>>(
        true,
        [&]() { return is_sorted_impl(an->lc, less); },
        [&]() { return is_sorted_impl(an->rc, less); });
    auto [l_sorted, l_left, l_right] = P.first;
    auto [r_sorted, r_left, r_right] = P.second;

    const ET& mid = Tree::get_entry(an);
    bool mid_sorted = less(l_right, mid) && less(mid, r_left);

    return {l_sorted && r_sorted && mid_sorted, l_left, r_right};
  }

  template <class F>
  static bool is_sorted(node* a, const F& less) {
    return std::get<0>(is_sorted_impl(a, less));
  }

  static std::optional<size_t> find_unsorted_compressed(node* a, const ET& e) {
    size_t i=0;
    std::optional<size_t> ret;
    auto fn = [&] (const ET& a) -> bool {
      if (e == a) {
        ret = i;
        return false;
      }
      i++;
      return true;
    };
    Tree::iterate_cond(a, fn);
    return ret;
  }

  static std::optional<size_t> find_unsorted(node* a, const ET& e, bool par=false) {
    if (a == nullptr) return {};
    if (Tree::is_compressed(a)) return find_unsorted_compressed(a, e);
    auto an = (regular_node*)a;

    auto P = utils::fork<std::optional<size_t>>(
        par,
        [&]() { return find_unsorted(an->lc, e, par); },
        [&]() { return find_unsorted(an->rc, e, true); });
    if (P.first.has_value()) return P.first;
    return P.second;
  }

  static node* reverse(node* a) {
    if (a == nullptr) return nullptr;
    if (Tree::is_compressed(a)) {
      ET arr[2*B];
      size_t k = Tree::size(a) - 1;
      auto fn = [&] (const ET& e) {
        arr[k] = e;
        k--;
      };
      Tree::iterate_seq(a, fn);
      return Tree::make_single_compressed_node(arr, Tree::size(a));
//      ET* arr = Tree::compressed_node_elms(a, tmp_arr);
//      size_t arr_size = Tree::size(a);
//      std::reverse(arr, arr + arr_size);
//      return to_tree_impl(arr, arr_size);
    }

    auto an = (regular_node*)a;
    size_t n = Tree::size(an);
    auto P = utils::fork<node*>(
        n >= kNodeLimit,
        [&]() { return reverse(an->lc); },
        [&]() { return reverse(an->rc); });

    const ET& mid = Tree::get_entry(an);
    regular_node* m = Tree::make_regular_node(mid);
    return Tree::node_join(P.second, P.first, m);
  }

  template <typename F>
  static void foreach_index(ptr a, size_t start, const F& f,
                            size_t granularity = kNodeLimit) {
    if (a.empty()) return;
    if (a.is_compressed()) {
      auto c = a.unsafe_ptr();  // Why unsafe ptr here?? node_ptr() causes ref_cnt bumps.
      size_t i = 0;
      auto fn = [&] (ET a) {
        f(a, start + i);
        i++;
      };
      Tree::iterate_seq(c, fn);
      return;
    }
    auto[lc, e, rc, root] = expose(std::move(a));
    size_t lsize = lc.size();
    f(e, start + lsize);
    utils::fork_no_result(
        lsize >= granularity,
        [&]() { foreach_index(std::move(lc), start, f, granularity); },
        [&]() {
          foreach_index(std::move(rc), start + lsize + 1, f, granularity);
        });
    GC::decrement(root);
  }

  // Experimental version: remove later if unused
  template <typename F>
  static void foreach_index_2(node* a, const F& f,
                            size_t granularity = kNodeLimit) {
    if (!a) return;
    if (Tree::is_compressed(a)) {
      auto fn = [&] (ET a) {
        f(a, 0);
      };
      Tree::iterate_seq(a, fn);
      return;
    }
    auto an = (regular_node*)a;
    f(Tree::get_entry(an), 0);
    utils::fork_no_result(true,
        [&]() { foreach_index_2(an->lc, f); },
        [&]() { foreach_index_2(an->rc, f); });
  }

  // F : entry x index -> bool
  template <typename F>
  static bool foreach_cond(node* a, const F& f) {
    if (a == nullptr) return true;
    if (Tree::is_compressed(a)) {
      return Tree::iterate_cond(a, f);
    }

    auto an = (regular_node*)a;
    if (foreach_cond(an->lc, f) && f(Tree::get_entry(a))) {
      return foreach_cond(an->rc, f);
    }
    return false;
  }

  // F : entry x index -> bool
  template <typename F, typename C>
  static bool foreach_cond_par(node* a, const F& f, const C& cond, bool par = false) {
    if (a == nullptr) return true;
    if (Tree::is_compressed(a)) {
      auto fn = [&] (const ET& a) -> bool {
        bool ret = f(a);
        return ret;
      };
      if (cond())
        return Tree::iterate_cond(a, fn);
      return false;
    }
    auto an = (regular_node*)a;
    bool ret = cond();
    if (ret) {
      ret = f(Tree::get_entry(an));
    }
    if (!ret) {
      return false;
    }
    auto P = utils::fork<bool>(
        par,
        [&]() { return foreach_cond_par(an->lc, f, cond, par); },
        [&]() {
          return foreach_cond_par(an->rc, f, cond, true);
    });
    return P.first && P.second;
  }

  // similar to above but sequential using in-order traversal
  // usefull if using 20th century constructs such as iterators
  template<typename F>
  static void foreach_seq(ptr a, const F& f) {
    if (a.empty()) return;
    auto[lc, e, rc, root] = expose(std::move(a));
    foreach_seq(std::move(lc), f);
    f(e);
    foreach_seq(std::move(rc), f);
    GC::decrement(root);
  }

  template <class Func>
  static node* filter_bc(ptr b1, const Func& f) {
    assert(b1.size() > 0);
    ET stack[kBaseCaseSize + 1];

    auto b1_node = b1.node_ptr();
    size_t offset = 0;
    auto copy_f = [&] (ET a) {  // needs to be a copy since we move
      if (f(a))
        parlay::move_uninitialized(stack[offset++], a);
    };
    Tree::iterate_seq(b1_node, copy_f);
    assert(offset <= kBaseCaseSize);

    Tree::decrement_recursive(b1_node);

    if (offset < B) {
      return to_tree_impl((ET*)stack, offset);
    } else {
      // need to refactor
      return Tree::make_compressed(stack, offset);
    }
  }

  template<class Func>
  static node* filter(ptr b, const Func& f, size_t granularity=kNodeLimit) {
    if (b.empty()) return NULL;
    size_t n = b.size();

//    if (n <= kBaseCaseSize) {
//      return filter_bc(std::move(b), f);
//    }

    if (b.is_compressed()) {
      return filter_bc(std::move(b), f);
    }

    auto[lc, e, rc, root] = expose(std::move(b));

    auto [l, r] = utils::fork<node*>(n >= granularity,
      [&]() {return filter(std::move(lc), f, granularity);},
      [&]() {return filter(std::move(rc), f, granularity);});

    if (f(e)) {
      if (!root) root = Tree::single(e);
      return join(l, e, r, root);
    } else {
      GC::decrement(root);
      return join2(l, r);
    }
  }

// TODO
//  template<class Func>
//  static bool if_exist(node* b, const Func& f, bool* indicator) {
//    if (!b) return false;
//	if (*indicator) return true;
//	if (f(Tree::get_entry(b))) {
//		utils::atomic_compare_and_swap(indicator, false, true);
//		//*indicator = true;
//		return true;
//	}
//
//    auto P = utils::fork<bool>(Tree::size(b) >= kNodeLimit,
//      [&]() {return if_exist(b->lc, f, indicator);},
//      [&]() {return if_exist(b->rc, f, indicator);});
//	return P.first || P.second;
//  }

  // Assumes the input is sorted and there are no duplicate keys
  static node* from_array(ET* A, size_t n) {
    if (n <= 0) return Tree::empty();
    if (n == 1) return Tree::single(A[0]);
    if (n >= B && n <= 2*B) {
      return Tree::make_compressed(A, n);
    }
    size_t mid = n/2;

    regular_node* m = Tree::make_regular_node(A[mid]);

    auto P = utils::fork<node*>(n >= kNodeLimit,
      [&]() {return from_array(A, mid);},
      [&]() {return from_array(A+mid+1, n-mid-1);});

    return Tree::node_join(P.first, P.second, m);
  }

// TODO
//  template<class Seq1, class Func>
//  static node* map_filter(typename Seq1::node* b, const Func& f,
//			  size_t granularity=kNodeLimit) {
//    if (!b) return NULL;
//
//    auto P = utils::fork<node*>(Seq1::size(b) >= granularity,
//      [&]() {return map_filter<Seq1>(b->lc, f, granularity);},
//      [&]() {return map_filter<Seq1>(b->rc, f, granularity);});
//
//    std::optional<ET> me = f(Seq1::get_entry(b));
//    if (me.has_value()) {
//      node* r = Tree::make_node(*me);
//      return Tree::node_join(P.first, P.second, r);
//    } else return join2(P.first, P.second);
//  }

  template<class R, class F>
  static typename R::T map_reduce(node* a, F f, R r,
				  size_t grain=kNodeLimit) {
    using T = typename R::T;
    if (a == nullptr) return r.identity();
    if (Tree::is_compressed(a)) {
      T v = r.identity();
      auto fn = [&] (ET a) {
        v = R::add(v, f(a));
      };
      Tree::iterate_seq(a, fn);
      return v;
    }

    size_t size = Tree::size(a);
    auto an = (regular_node*)a;
    auto P = utils::fork<T>(size >= grain,
      [&]() {return map_reduce<R>(an->lc, f, r, grain);},
      [&]() {return map_reduce<R>(an->rc, f, r, grain);});

    T v = f(Tree::get_entry(a));
    return R::add(P.first, r.add(v, P.second));
  }

// TODO
//  template<class F, class T>
//  static void semi_map_reduce_seq(node* b, T& v, F& f) {
//    if (!b) return;
//    semi_map_reduce_seq(b->lc, v, f);
//    f(v, Tree::get_entry(b));
//    semi_map_reduce_seq(b->rc, v, f);
//  }

// TODO
//  // Class R must be a structure representing a monoid with:
//  //   T = a type
//  //   T identity();
//  //   T add(T, T);  an associative operation on T
//  // Function f must be of type (T&, E) -> void
//  //   it side affects its first argument by adding to it
//  template<class R, class F>
//  static typename R::T semi_map_reduce(node* b, F& f, R reduce, size_t grain) {
//    using T = typename R::T;
//    if (Tree::size(b) >= grain) {
//      T r, l;
//      auto left = [&] () -> void {
//	r = semi_map_reduce<R,F>(b->rc, f, reduce, grain);};
//      auto right = [&] () -> void {
//	l = semi_map_reduce<R,F>(b->lc, f, reduce, grain);};
//      parlay::par_do(left,right);
//      f(l, Tree::get_entry(b));
//      return reduce.add(l,r);
//    } else {
//      // when going sequential, create identity and then update it
//      T id = reduce.identity();
//      semi_map_reduce_seq(b, id, f);
//      return id;
//    }
//  }

};

}  // namespace cpam
