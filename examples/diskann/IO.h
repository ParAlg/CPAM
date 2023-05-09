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

// *************************************************************
//  BINARY TOOLS: uint8, int8, float32, int32
// *************************************************************

auto parse_uint8bin(const char* filename){
    auto [fileptr, length] = mmapStringFromFile(filename);

    int num_vectors = *((int*) fileptr);
    int d = *((int*) (fileptr+4));

    std::cout << "Detected " << num_vectors << " points with dimension " << d << std::endl;
    parlay::sequence<Tvec_point<uint8_t>> points(num_vectors);

    parlay::parallel_for(0, num_vectors, [&] (size_t i) {
        points[i].id = i; 

        uint8_t* start = (uint8_t*)(fileptr + 8 + i*d); //8 bytes at the start for size + dimension
        uint8_t* end = start + d;
        points[i].coordinates = parlay::make_slice(start, end);
    });

    return points;
}

auto parse_int8bin(const char* filename){
    auto [fileptr, length] = mmapStringFromFile(filename);

    int num_vectors = *((int*) fileptr);
    int d = *((int*) (fileptr+4));

    std::cout << "Detected " << num_vectors << " points with dimension " << d << std::endl;
    parlay::sequence<Tvec_point<int8_t>> points(num_vectors);

    parlay::parallel_for(0, num_vectors, [&] (size_t i) {
        points[i].id = i; 

        int8_t* start = (int8_t*)(fileptr + 8 + i*d); //8 bytes at the start for size + dimension
        int8_t* end = start + d;
        points[i].coordinates = parlay::make_slice(start, end);
    });

    return points;
}

auto parse_fbin(const char* filename){
    auto [fileptr, length] = mmapStringFromFile(filename);

    int num_vectors = *((int*) fileptr);
    int d = *((int*) (fileptr+4));

    std::cout << "Detected " << num_vectors << " points with dimension " << d << std::endl;
    parlay::sequence<Tvec_point<float>> points(num_vectors);

    parlay::parallel_for(0, num_vectors, [&] (size_t i) {
        points[i].id = i; 

        float* start = (float*)(fileptr + 8 + 4*i*d); //8 bytes at the start for size + dimension
        float* end = start + d;
        points[i].coordinates = parlay::make_slice(start, end);
    });

    return points;
}

auto parse_ibin(const char* filename){
    auto [fileptr, length] = mmapStringFromFile(filename);

    int num_vectors = *((int*) fileptr);
    int d = *((int*) (fileptr+4));

    std::cout << "Detected " << num_vectors << " points with number of results " << d << std::endl;
    parlay::sequence<ivec_point> points(num_vectors);

    parlay::parallel_for(0, num_vectors, [&] (size_t i) {
        points[i].id = i; 

        int* start = (int*)(fileptr + 8 + 4*i*d); //8 bytes at the start for size + dimension
        int* end = start + d;
        float* dist_start = (float*)(fileptr+ 8 + num_vectors*4*d + 4*i*d);
        float* dist_end = dist_start+d; 
        points[i].coordinates = parlay::make_slice(start, end);
        points[i].distances = parlay::make_slice(dist_start, dist_end);
    });

    return points;
}
