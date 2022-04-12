#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <parallel/algorithm>

#include <cpam/cpam.h>
#include <pam/pam.h>
#include <pam/get_time.h>

timer init_tm, build_tm, total_tm, aug_tm, sort_tm, reserve_tm, outmost_tm, globle_tm, g_tm;
double early;

#ifdef USE_PAM
#define ENTRY_TYPE pair
#else
#define ENTRY_TYPE tuple
#endif

template<typename c_type, typename w_type>
struct RangeQuery {
  using x_type = c_type;
  using y_type = c_type;
  using point_type = Point<c_type, w_type>;
  using point_x = pair<x_type, y_type>;
  using point_y = pair<y_type, x_type>;

  struct inner_map_t {
    using key_t = point_y;
    using val_t = w_type;
    static bool comp(key_t a, key_t b) { return a < b;}
    using aug_t = w_type;
    inline static aug_t from_entry(key_t k, val_t v) { return v; }
    inline static aug_t combine(aug_t a, aug_t b) { return a+b; }
    static aug_t get_empty() { return 0;}
  };
#ifdef USE_PAM
  using inner_map = aug_map<inner_map_t>;
#else
  using inner_map = cpam::aug_map<inner_map_t, 32>;
#endif

  struct outer_map_t {
    using key_t = point_x;
    using val_t = w_type;
    static bool comp(key_t a, key_t b) { return a < b;}
    using aug_t = inner_map;
    inline static aug_t from_entry(key_t k, val_t v) {
      return aug_t({ENTRY_TYPE(point_y(k.second, k.first), v)});
    }
    inline static aug_t combine(aug_t a, aug_t b) {
      auto ret = aug_t::map_union(std::move(a), std::move(b), [](w_type x, w_type y) {return x+y;});
      return ret;
    }
    static aug_t get_empty() { return aug_t();}
  };
#ifdef USE_PAM
  using outer_map = aug_map<outer_map_t>;
#else
  using outer_map = cpam::aug_map<outer_map_t, 128>;
#endif

  RangeQuery(sequence<point_type>& points) {
    std::cout << "Calling construct." << std::endl;
    construct(points);
    std::cout << "Finished construct. " << std::endl;
  }

  ~RangeQuery() {
  }

  static void reserve(size_t n) {
    reserve_tm.start();
    outer_map::reserve(n/64);
    inner_map::reserve(24*n/64);
    reserve_tm.stop();
  }

  static void finish() {
    outer_map::finish();
    inner_map::finish();
  }

  void construct(sequence<point_type>& points) {
    const size_t n = points.size();
    reserve(n);
    total_tm.start();
#ifdef USE_PAM
    using my_pair = std::pair<point_x, w_type>;
#else
    using my_pair = std::tuple<point_x, w_type>;
#endif

    auto pointsEle = tabulate(n, [&] (size_t i) -> my_pair {
      return my_pair(make_pair(points[i].x, points[i].y), points[i].w);
    });
    std::cout << "Built pointsEle" << std::endl;
    range_tree = outer_map(pointsEle);

    total_tm.stop();
  }


  struct count_t {
    point_y y1, y2;
    int r;
    count_t(point_y y1, point_y y2) : y1(y1), y2(y2), r(0) {}
    void add_entry(ENTRY_TYPE<point_x,w_type> e) {
      if (std::get<0>(e).second >= y1.first && std::get<0>(e).second <= y2.first) r += 1;
    }
    void add_aug_val(inner_map a) {
      auto inc = (a.rank(y2) - a.rank(y1));
      r += inc;}
  };


  w_type query_count(x_type x1, y_type y1, x_type x2, y_type y2) {
    count_t qrc(make_pair(y1,x1), make_pair(y2,x2));
    range_tree.range_sum(make_pair(x1,y1), make_pair(x2,y2), qrc);
    return qrc.r;
  }

  void check_all() {
    size_t ct = 0;
    auto iter_f = [&] (const auto& et, const auto& aug) {
      point_x key = et.first;
      auto& inner = et.second;
      std::cout << "Entry " << ct << " key = (" << key.first << "," << key.second << ") " << inner << " aug_size = " << aug.size() << " aug_ref = " << aug.ref_cnt() << std::endl;
      ct++;
    };
    range_tree.iterate_aug(iter_f);
  }

  struct sum_t {
    point_y y1, y2;
    int r;
    sum_t(point_y y1, point_y y2) : y1(y1), y2(y2), r(0) {}
    void add_entry(ENTRY_TYPE<point_x,w_type> e) {
      if (std::get<0>(e).second >= y1.first && std::get<0>(e).second <= y2.first) r += std::get<1>(e);
    }
    void add_aug_val(inner_map a) { r += a.aug_range(y1, y2); }
  };

  w_type query_sum(x_type x1, y_type y1, x_type x2, y_type y2) {
    sum_t qrs(make_pair(y1,x1), make_pair(y2,x2));
    range_tree.range_sum(make_pair(x1,y1), make_pair(x2,y2), qrs);
    return qrs.r;
  }

  struct range_t {
    point_y y1, y2;
    point_y* out;
    range_t(point_y y1, point_y y2, point_y* out) : y1(y1), y2(y2), out(out) {}
    void add_entry(ENTRY_TYPE<point_x,w_type> e) {
      if (std::get<0>(e).second >= y1.first && std::get<0>(e).second <= y2.first) {
        *out++ = point_y(std::get<0>(e).second,std::get<0>(e).second);
      }
    }
    void add_aug_val(inner_map a) {
      inner_map r = inner_map::range(a, y1, y2);
      size_t s = r.size();
      inner_map::keys(r,out);
      out += s;
    }
  };

  sequence<point_y> query_range(x_type x1, y_type y1, x_type x2, y_type y2) {
    // std::cout << "Calling query count on " << x1 << " " << y1 << " " << x2 << " " << y2 << std::endl;
    size_t n = query_count(x1,y1,x2,y2);
    // std::cout << "Count result = " << n << std::endl;
    auto out = sequence<point_y>::uninitialized(n);
    range_t qr(point_y(y1, x1), point_y(y2, x2), out.data());

    // std::cout << "Calling range_sum" << std::endl;
    range_tree.range_sum(point_x(x1, y1), point_x(x2, y2), qr);
    // std::cout << "Finished range_sum" << std::endl;
    return out;
  }

  void insert_point(point_type p) {
    range_tree.insert(make_pair(make_pair(p.x, p.y), p.w));
  }

  void insert_point_lazy(point_type p) {
    range_tree = outer_map::insert_lazy(std::move(range_tree),
					make_pair(make_pair(p.x, p.y), p.w));
  }

  void print_status() {
    cout << "Outer Map: ";  outer_map::GC::print_stats();
    cout << "Inner Map: ";  inner_map::GC::print_stats();
    parlay::internal::get_default_allocator().print_stats();
  }

  size_t size_in_bytes() {
    // PAM and CPAM-internal specific nodes:
    size_t outer_used = outer_map::GC::used_node();
    size_t inner_used = inner_map::GC::used_node();
#ifdef USE_PAM
    size_t internal_nodes_space = (sizeof(typename outer_map::node) * outer_used) + (sizeof(typename inner_map::node) * inner_used);
#else
    size_t internal_nodes_space = (sizeof(typename outer_map::GC::regular_node) * outer_used) + (sizeof(typename inner_map::GC::regular_node) * inner_used);
#endif
    auto [used, unused] = parlay::internal::get_default_allocator().stats();
    return used + internal_nodes_space;
  }

  bool find(point_type p) {
    return range_tree.find(make_pair(p.x, p.y));
  }

  outer_map range_tree;
};
