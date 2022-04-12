#pragma once

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <stdio.h>
#include <string.h>

#include "macros.h"
#include "compression.h"

namespace aspen {
namespace bytepd_amortized {

#define PARALLEL_DEGREE 1000

  inline size_t get_virtual_degree(uintE d, uchar* ngh_arr) {
    if (d > 0) {
      return *((uintE*)ngh_arr);
    }
    return 0;
  }

  inline uintE eatFirstEdge(uchar* &start, uintE source) {
    uchar fb = *start++;
    uintE edgeRead = (fb & 0x3f);
    if (LAST_BIT_SET(fb)) {
      int shiftAmount = 6;
      while (1) {
        uchar b = *start;
        edgeRead |= ((b & 0x7f) << shiftAmount);
        start++;
        if (LAST_BIT_SET(b))
          shiftAmount += EDGE_SIZE_PER_BYTE;
        else
          break;
      }
    }
    return (fb & 0x40) ? source - edgeRead : source + edgeRead;
  }

  inline uintE eatEdge(uchar* &start) {
    uintE edgeRead = 0;
    int shiftAmount = 0;

    while (1) {
      uchar b = *start;
      edgeRead += ((b & 0x7f) << shiftAmount);
      start++;
      if (LAST_BIT_SET(b))
        shiftAmount += EDGE_SIZE_PER_BYTE;
      else
        break;
    }
    return edgeRead;
  }

  struct simple_iter {
    uchar* base;
    uchar* finger;
    uintE src;
    uintE degree;

    uintE num_blocks;
    uintE cur_chunk;
    uintE last_edge;

    uintE proc;

    simple_iter(uchar* _base, uintE _degree, uintE _src) : base(_base), src(_src), degree(_degree), 
                                                    cur_chunk(0) {
      if (degree == 0) return;
      uintE virtual_degree = *((uintE*)base);
      num_blocks = 1+(virtual_degree-1)/PARALLEL_DEGREE;

      finger = base + (num_blocks-1)*sizeof(uintE) + sizeof(uintE);

      finger += sizeof(uintE);
      last_edge = eatFirstEdge(finger, src);
      proc = 1;
    }

    inline uintE cur() {
      return last_edge;
    }

    inline uintE next() {
      if (proc == PARALLEL_DEGREE) {
        finger += sizeof(uintE); // skip block start
        last_edge = eatFirstEdge(finger, src);
        proc = 1;
        cur_chunk++;
      } else {
        last_edge += eatEdge(finger);
        proc++;
      }
      return last_edge;
    }

    inline bool has_next() {
      return (cur_chunk*PARALLEL_DEGREE + proc) < degree;
    }
  };

}  // namespace bytepd_amortized
}  // namespace aspen
