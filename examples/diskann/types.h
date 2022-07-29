#pragma once

#include <algorithm>
#include "parlay/parallel.h"
#include "parlay/primitives.h"

using node_id = uint32_t;

// for a file in .fvecs or .bvecs format, but extendible to other types
template <typename T>
struct Tvec_point {
  node_id id;
  parlay::slice<T*, T*> coordinates;
  Tvec_point() : coordinates(parlay::make_slice<T*, T*>(nullptr, nullptr)) {}
};

// for an ivec file, which contains the ground truth
// only info needed is the coordinates of the nearest neighbors of each point
struct ivec_point {
  node_id id;
  parlay::slice<int*, int*> coordinates;
  ivec_point()
      : coordinates(parlay::make_slice<int*, int*>(nullptr, nullptr)) {}
};
