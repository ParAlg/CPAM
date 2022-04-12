#pragma once

#include <limits.h>
#include "aspen/aspen.h"

namespace aspen {

// a 32-bit hash function
inline uint32_t hash32(uint32_t a) {
  a = (a+0x7ed55d16) + (a<<12);
  a = (a^0xc761c23c) ^ (a>>19);
  a = (a+0x165667b1) + (a<<5);
  a = (a+0xd3a2646c) ^ (a<<9);
  a = (a+0xfd7046c5) + (a<<3);
  a = (a^0xb55a4f09) ^ (a>>16);
  return a;
}

template <class intT>
struct rMat {
  double a, ab, abc;
  intT n;
  intT h;
  rMat(intT _n, intT _seed,
       double _a, double _b, double _c) {
    n = _n; a = _a; ab = _a + _b; abc = _a+_b+_c;
    h = hash32((intT)_seed);
    if(abc > 1) { cout << "in rMat: a + b + c add to more than 1\n"; abort();}
    if((1UL << parlay::log2_up(n)) != n) { cout << "in rMat: n not a power of 2"; abort(); }
  }


  double hashDouble(intT i) {
    return ((double) (hash32((intT)i))/((double) std::numeric_limits<int32_t>::max()));}

  std::pair<intT, intT> rMatRec(intT nn, intT randStart, intT randStride) {
    if (nn==1) return std::make_pair<intT, intT>(0,0);
    else {
      std::pair<intT, intT> x = rMatRec(nn/2, randStart + randStride, randStride);
      double r = hashDouble(randStart);
      if (r < a) return x;
      else if (r < ab) return std::make_pair(x.first,x.second+nn/2);
      else if (r < abc) return std::make_pair(x.first+nn/2, x.second);
      else return std::make_pair(x.first+nn/2, x.second+nn/2);
    }
  }

  std::pair<intT, intT> operator() (intT i) {
    intT randStart = hash32((intT)(2*i)*h);
    intT randStride = hash32((intT)(2*i+1)*h);
    return rMatRec(n, randStart, randStride);
  }
};

}  // namespace apsne
