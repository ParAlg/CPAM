#pragma once

#include <iostream>
#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "cpam/parse_command_line.h"
#include "types.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// returns a pointer and a length
std::pair<char*, size_t> mmapStringFromFile(const char* filename) {
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

auto parse_fvecs(const char* filename, int maxDeg) {
  auto [fileptr, length] = mmapStringFromFile(filename);

  // Each vector is 4 + 4*d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are floats representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);

  size_t vector_size = 4 + 4*d;
  size_t num_vectors = length / vector_size;
  // std::cout << "Num vectors = " << num_vectors << std::endl;

  parlay::sequence<Tvec_point<float>> points(num_vectors);

  parlay::sequence<int> &out_nbh = *new parlay::sequence<int>(maxDeg*num_vectors);

  parlay::parallel_for(0, num_vectors*maxDeg, [&] (size_t i){out_nbh[i] = -1; });

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    float* start = (float*)(fileptr + offset_in_bytes);
    float* end = start + d;
    points[i].id = i;
    points[i].coordinates = parlay::make_slice(start, end);
    points[i].out_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));
    // points[i].new_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));
  });

  return points;
}

auto parse_ivecs(const char* filename) {
  auto [fileptr, length] = mmapStringFromFile(filename);

  // Each vector is 4 + 4*d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are floats representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);

  size_t vector_size = 4 + 4*d;
  size_t num_vectors = length / vector_size;

  parlay::sequence<ivec_point> points(num_vectors);

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    int* start = (int*)(fileptr + offset_in_bytes);
    int* end = start + d;
    points[i].id = i;
    points[i].coordinates = parlay::make_slice(start, end);
  });

  return points;
}

auto parse_bvecs(const char* filename, int maxDeg) {

  auto [fileptr, length] = mmapStringFromFile(filename);
  // Each vector is 4 + d bytes.
  // * first 4 bytes encode the dimension (as an integer)
  // * next d values are unsigned chars representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  int d = *((int*)fileptr);
  // std::cout << "Dimension = " << d << std::endl;

  size_t vector_size = 4 + d;
  size_t num_vectors = length / vector_size;

  parlay::sequence<Tvec_point<uint8_t>> points(num_vectors);

  parlay::sequence<int> &out_nbh = *new parlay::sequence<int>(maxDeg*num_vectors);

  parlay::parallel_for(0, num_vectors*maxDeg, [&] (size_t i){out_nbh[i] = -1;});

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    size_t offset_in_bytes = vector_size * i + 4;  // skip dimension
    uint8_t* start = (uint8_t*)(fileptr + offset_in_bytes);
    uint8_t* end = start + d;
    points[i].id = i;
    points[i].coordinates = parlay::make_slice(start, end);
    points[i].out_nbh = parlay::make_slice(out_nbh.begin()+maxDeg*i, out_nbh.begin()+maxDeg*(i+1));
  });

  return points;
}

