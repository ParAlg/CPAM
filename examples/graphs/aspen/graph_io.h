#pragma once

#include "macros.h"

#include "parlay/io.h"

#include <fcntl.h>
#if defined(__APPLE__)
#else
#include <malloc.h>
#endif
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <cmath>
#include <fstream>

#include "compression.h"

namespace aspen {

template <class weight_type>
struct Edge {
  uintE from;
  uintE to;
  weight_type weight;

  Edge() {}
  Edge(uintE _from, uintE _to)
    : from(_from)
    , to(_to)
    , weight(weight_type()) {}
  Edge(const uintE _from, const uintE _to, const weight_type _weight)
    : from(_from)
    , to(_to)
    , weight(_weight) {}
};

namespace internal {  // Internal declarations

// Header string expected at the top of unweighted adjacency graph files.
const std::string kUnweightedAdjGraphHeader = "AdjacencyGraph";
// Header string expected at the top of weighted adjacency graph files.
const std::string kWeightedAdjGraphHeader = "WeightedAdjacencyGraph";

} // namespace internal

// returns a pointer and a length
inline std::pair<char*, size_t> mmapStringFromFile(const char* filename) {
  struct stat sb;
  int fd = open(filename, O_RDONLY);
  if (fd == -1) {
    perror("open");
    exit(-1);
  }
  if (fstat(fd, &sb) == -1) {
    perror("fstat");
    exit(-1);
  }
  if (!S_ISREG(sb.st_mode)) {
    perror("not a file\n");
    exit(-1);
  }
  char* p =
      static_cast<char*>(mmap(0, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0));
  if (p == MAP_FAILED) {
    perror("mmap");
    exit(-1);
  }
  if (close(fd) == -1) {
    perror("close");
    exit(-1);
  }
  size_t n = sb.st_size;
  return std::make_pair(p, n);
}

inline parlay::sequence<char> readStringFromFile(const char* fileName) {
  std::ifstream file(fileName, std::ios::in | std::ios::binary | std::ios::ate);
  if (!file.is_open()) {
    abort();
  }
  uint64_t end = file.tellg();
  file.seekg(0, std::ios::beg);
  uint64_t n = end - file.tellg();
  auto bytes = parlay::sequence<char>(n); // n+1?
  file.read(bytes.begin(), n);
  file.close();
  return bytes;
}


/* Returns a tuple containing (n, m, offsets, edges) --- the number of
 * vertices, edges, the vertex offsets, and the edge values, after
 * parsing the input graph file */
inline std::tuple<size_t, size_t, uintT*, uintE*> parse_unweighted_symmetric_graph(
    const char* fname,
    bool mmap,
    char* bytes = nullptr,
    size_t bytes_size = 0) {
  // parlay::sequence<parlay::slice<char*, char*>> tokens;
  parlay::sequence<char> S;

  if (bytes == nullptr) {
    if (mmap) {
      std::pair<char*, size_t> MM = mmapStringFromFile(fname);
      S = parlay::sequence<char>::uninitialized(MM.second);
      parlay::parallel_for(0, S.size(), [&] (size_t i) { S[i] = MM.first[i]; });
      if (munmap(MM.first, MM.second) == -1) {
        perror("munmap");
        exit(-1);
      }
    } else {
      S = readStringFromFile(fname);
    }
  }
  parlay::sequence<parlay::slice<char*, char*>> tokens = parlay::map_tokens(parlay::make_slice(S),[] (auto x) { return parlay::make_slice(x); });

  uint64_t n = parlay::internal::chars_to_int_t<unsigned long>(tokens[1]);
  uint64_t m = parlay::internal::chars_to_int_t<unsigned long>(tokens[2]);

  uintT* offsets = cpam::utils::new_array_no_init<uintT>(n+1);
  uintE* edges = cpam::utils::new_array_no_init<uintE>(m);

  parlay::parallel_for(0, n, [&] (size_t i)
                  { offsets[i] = parlay::internal::chars_to_int_t<unsigned long>(tokens[i + 3]); });
  offsets[n] = m; /* make sure to set the last offset */
  parlay::parallel_for(0, m, [&] (size_t i)
                  { edges[i] = parlay::internal::chars_to_int_t<unsigned long>(tokens[i + n + 3]); });
  S.clear();

  tokens.clear();
  return std::make_tuple(n, m, offsets, edges);
}

auto read_o_direct(const char* fname) {
  int fd;
  if ( (fd = open(fname, O_RDONLY | O_DIRECT) ) != -1) {
    cout << "input opened!" << endl;
  } else {
    cout << "can't open input file!";
  }

  size_t fsize = lseek(fd, 0, SEEK_END);
  lseek(fd, 0, 0);
  auto s = (char*)memalign(4096 * 2, fsize + 4096);

  cout << "fsize = " << fsize << endl;

  size_t sz = 0;

  size_t pgsize = getpagesize();
  cout << "pgsize = " << pgsize << endl;

  size_t read_size = 1024*1024*1024;
  if (sz + read_size > fsize) {
    size_t k = std::ceil((fsize - sz) / pgsize);
    read_size = std::max(k*pgsize, pgsize);
    cout << "set read size to: " << read_size << " " << (fsize - sz) << " bytes left" << endl;
  }

  while (sz + read_size < fsize) {
    void* buf = s + sz;
    cout << "reading: " << read_size << endl;
    sz += read(fd, buf, read_size);
    cout << "read: " << sz << " bytes" << endl;
    if (sz + read_size > fsize) {
      size_t k = std::ceil((fsize - sz) / pgsize);
      read_size = std::max(k*pgsize, pgsize);
      cout << "set read size to: " << read_size << " " << (fsize - sz) << " bytes left" << endl;
    }
  }
  if (sz < fsize) {
    cout << "last read: rem = " << (fsize - sz) << endl;
    void* buf = s + sz;
    sz += read(fd, buf, pgsize);
    cout << "read " << sz << " bytes " << endl;
  }
  return s;
}

auto parse_unweighted_compressed_symmetric_graph(const char* fname, bool mmap=false) {
  char* s;
  if (mmap) {
    auto SS = mmapStringFromFile(fname);
    s = SS.first;
  } else {
    s = read_o_direct(fname);
  }

  return s;
}


}  // namespace aspen
