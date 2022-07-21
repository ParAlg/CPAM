#pragma once

#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"

// for a file in .fvecs or .bvecs format, but extendible to other types
template <typename T>
struct Tvec_point {
  int id;
  int cnt;
  parlay::slice<T*, T*> coordinates;
  parlay::slice<int*, int*> out_nbh;
  parlay::slice<int*, int*> new_nbh;
  Tvec_point()
      : coordinates(parlay::make_slice<T*, T*>(nullptr, nullptr)),
        out_nbh(parlay::make_slice<int*, int*>(nullptr, nullptr)),
        new_nbh(parlay::make_slice<int*, int*>(nullptr, nullptr)) {}
  parlay::sequence<int> ngh = parlay::sequence<int>();
};

// for an ivec file, which contains the ground truth
// only info needed is the coordinates of the nearest neighbors of each point
struct ivec_point {
  int id;
  parlay::slice<int*, int*> coordinates;
  ivec_point()
      : coordinates(parlay::make_slice<int*, int*>(nullptr, nullptr)) {}
};
