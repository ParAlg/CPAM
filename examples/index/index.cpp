#include <algorithm>
#include <iostream>
#include <cstring>
#include <vector>
#include <cctype>
#include <fstream>
#include <set>

#include <parlay/io.h>
#include <parlay/primitives.h>
#include <parlay/random.h>

#include "index.h"

using namespace std;

long str_to_long(char* str) {
    return strtol(str, NULL, 10);
}

using post_elt = inv_index::post_elt;
using post_list = inv_index::post_list;
using index_elt = inv_index::index_elt;
using index_pair = std::tuple<token,post_list>;

template <class Seq>
parlay::sequence<index_elt> parse(Seq const &Str, bool verbose, timer t) {
  size_t n = Str.size();

  // find start of each document based on header
  string header = "<doc id=";
  size_t header_size = header.size();
  auto doc_ends = parlay::tabulate(n, [&] (size_t i) -> bool {
      if (i > n - header_size - 1) return false;
      for (size_t j=0; j < header_size; j++)
	if (header[j] != Str[i+j+1]) return false;
      return true;
    });
  t.next("get headers");

  // only consider alphabetic characters, and turn to lowercase
  parlay::sequence<char> cleanstr = parlay::map(Str, [&] (char a) -> char {
      return isalpha(a) ? tolower(a) : ' ';});
  t.next("clean");

  // partition cleaned string based on document ends
  auto docs = parlay::split_at(cleanstr, doc_ends);
  std::cout << "Num docs = " << docs.size() << std::endl;
  t.next("partition");

  // tokenize each document and tag with document id and weight
  auto pairs = parlay::tabulate(docs.size(), [&] (size_t i) {
      //auto is_space = [] (char a) {return a == ' ' || a == '\t';};
      auto is_space = [] (char a) {return isspace(a);};
      size_t m = docs[i].size();
      // remove header from each doc
      auto doc = docs[i].cut(std::min(header_size, m), m);

      // get tokens from the remaining doc, and tag each with doc id and weight
      // not sure if map or dmap is faster
      return parlay::dmap(parlay::tokens(doc, is_space), [=] (parlay::sequence<char> w) {
	  return index_elt(w, post_elt(i,1.0));});
    });
  t.next("tokens");

  auto r = parlay::flatten(pairs);
  t.next("flatten");

  if (verbose) {
    cout << "Words = " << r.size() << endl;
    cout << "Documents = " << docs.size() << endl;
  }

  return r;
}

int main(int argc, char** argv) {
  string default_filename = "wiki_small.txt";
  commandLine P(argc, argv,
		"./index [-o] [-v] [-r rounds] [-n max_chars] [-q num_queries] [-f file]");
  size_t max_chars = P.getOptionLongValue("-n", 1000000000000);
  size_t num_queries = P.getOptionLongValue("-q", 10000);
  string fname = string(P.getOptionValue("-f", default_filename));
  bool verbose = P.getOption("-v");
  int rounds = P.getOptionIntValue("-r", 1);
  size_t threads = parlay::num_workers();
  timer t;
  timer tdetail("build index detail", verbose);

  // read
  auto Str_ = parlay::file_map(fname);
  auto Str__ = parlay::make_slice(Str_);
  tdetail.next("mmap file");
  auto Str = Str__.cut(0,std::min(max_chars,Str__.size()));
  tdetail.next("copy file");

  // parse input
  parlay::sequence<index_elt> KV = parse(Str, verbose, tdetail);
  size_t n = KV.size();
  inv_index test_idx;

  parlay::sequence<double> build_times;
  size_t index_size = 0;
  for (int i=0; i < rounds; i++) {
    test_idx = inv_index();
    // Build the index
    t.start();
    test_idx = inv_index(KV);
    double t_build = t.stop();
    build_times.push_back(t_build);
    if (index_size == 0) {
      index_size = test_idx.index_size();
    }
    if (verbose)
      cout << "unique words = " << test_idx.idx.size() << endl;
  }

  auto less = [] (double a, double b) {return a < b;};
  parlay::sort(parlay::make_slice(build_times), less);

  std::cout << "RESULT, name = index build, threads = " << threads
	 << ", rounds = 1"
	 << ", n = " << n
	 << ", q = " << num_queries
	 << ", time = " << build_times[rounds/2]
   << ", size_in_GiB = " << (static_cast<double>(index_size) / std::pow(1024, 3))
	 << endl;

  if (verbose) {
    parlay::internal::memory_usage();
    std::cout << "Per-bucket details: " << std::endl;
    parlay::internal::get_default_allocator().print_stats();
  }

  using idx = typename inv_index::index;

  parlay::sequence<pair<token,token>> test_word_pairs(num_queries);
  size_t* size_in = new size_t[num_queries];
  size_t* size_out = new size_t[num_queries];

  // Setup the queries
  // filters out any words which appear fewer than 100 times
  idx common_words = idx::filter(test_idx.idx, [] (index_pair e) {
      return (std::get<1>(e).size() > 100); });
  if (verbose)
    cout << "filter size = " << common_words.size() << endl;

  //size_t num_words = test_idx.idx.size();
  size_t num_common_words = common_words.size();
  if (num_common_words == 0) {
    cout << "data too small for running query" << endl;
    return 0;
  }

  parlay::random r(0);
  for(size_t i =0; i < num_queries; i++) {
    test_word_pairs[i].first = std::get<0>(KV[r.ith_rand(i) % (n - 1)]);
    auto rank = r.ith_rand(num_queries+i)%num_common_words;
    test_word_pairs[i].second =
      std::get<0>((*common_words.select(rank)));
  }

  // run the queries
  parlay::sequence<double> queries;
  for (int i=0; i < rounds; i++) {
    t.start();
    parlay::parallel_for(0, num_queries, [&] (size_t i) {
      post_list l1 = test_idx.get_list(test_word_pairs[i].first);
      post_list l2 = test_idx.get_list(test_word_pairs[i].second);
      size_in[i] = l1.size() + l2.size();
      post_list l3 = test_idx.And(std::move(l1),std::move(l2));
      vector<post_elt> r = test_idx.top_k(l3,10);
      size_out[i] = r.size();
    }, 1);
    double t_query = t.stop();
    queries.push_back(t_query);
  }
  size_t total_in = 0;
  size_t total_out = 0;

  parlay::sort(parlay::make_slice(queries), less);

  std::cout << "RESULT, name = index query, threads = " << threads
       << ", rounds = 1"
       << ", n = " << n
       << ", q = " << num_queries
       << ", time = " << queries[rounds/2]
       << endl;

  for (size_t i=0; i < num_queries; i++) {
    total_in += size_in[i];
    total_out += size_out[i];
  }

  if (verbose) {
    cout << "total in = " << total_in << endl;
    cout << "total out = " << total_out << endl;

    cout << "Unique words = " << test_idx.idx.size() << endl;
    cout << "Map: ";  inv_index::index::GC::print_stats();
    cout << "Set: ";  post_list::GC::print_stats();

    test_idx.print_index_size();
  }

  return 0;
}
