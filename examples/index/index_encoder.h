#pragma once

#include <functional>
#include <optional>
#include <tuple>
#include <type_traits>

#include "cpam/byte_encode.h"

namespace cpam {

struct diffencoded_index_encoder {

  struct data {};

  template <class Entry, bool is_aug = false>
  struct encoder {

    using ET = typename Entry::entry_t;
    using K = int;  // inv_index::doc_id
    using V = int;  // inv_index::weight

    static inline void print_info(const ET& et) {}

    static inline size_t encoded_size(ET* data, size_t size) {
      uint8_t stk[2*sizeof(K)];
      assert(size > 0);
      K prev_key = Entry::get_key(data[0]);
      size_t key_bytes = sizeof(K);  // first key is uncompressed
      size_t val_bytes = 0;
      for (size_t i=0; i<size; i++) {
        K cur_key = Entry::get_key(data[i]);
        K next_diff = cur_key - prev_key;
        auto bytes_used = encodeUnsigned<K>((uint8_t*)stk, 0, next_diff);  // compressed difference
        key_bytes += bytes_used;
        prev_key = cur_key;

        V cur_val = Entry::get_val(data[i]);
        val_bytes += encodeUnsigned((uint8_t*)stk, 0, cur_val);
      }
      return key_bytes + val_bytes;
    }

    static inline auto encode(ET* data, size_t size, uint8_t* bytes) {
      uint8_t* start = bytes;

      if constexpr (is_aug) {
        using AT = typename Entry::aug_t;
        AT av = Entry::from_entry(data[0]);
        K prev_key = Entry::get_key(data[0]);
        *((K*)start) = prev_key;  // store first key
        size_t offset = sizeof(K);  // first key is uncompressed
        offset = encodeUnsigned<V>(start, offset, Entry::get_val(data[0]));  // encode first val.
        for (size_t i=1; i<size; i++) {
          K cur_key = Entry::get_key(data[i]);
          K next_diff = cur_key - prev_key;
          offset = encodeUnsigned<K>(start, offset, next_diff);  // compressed difference
          offset = encodeUnsigned<V>(start, offset, Entry::get_val(data[i]));  // next val
          prev_key = cur_key;
          av = Entry::combine(av, Entry::from_entry(data[i]));
        }
        return av;
      } else {
        K prev_key = Entry::get_key(data[0]);
        *((K*)start) = prev_key;  // store first key
        size_t offset = sizeof(K);  // first key is uncompressed
        offset = encodeUnsigned<V>(start, offset, Entry::get_val(data[0]));  // encode first val.
        for (size_t i=1; i<size; i++) {
          K cur_key = Entry::get_key(data[i]);
          assert(cur_key > prev_key);
          K next_diff = cur_key - prev_key;
          offset = encodeUnsigned<K>(start, offset, next_diff);  // compressed difference
          offset = encodeUnsigned<V>(start, offset, Entry::get_val(data[i]));  // next val
          prev_key = cur_key;
        }
      }
    }

    template <class F>
    static inline void decode(uint8_t* bytes, size_t size, const F& f) {
      uint8_t* start = bytes;

      K prev_key = *((K*)start);
      start += sizeof(K);
      V cur_val = decodeUnsigned<V>(start);
      f(Entry::to_entry(prev_key, cur_val));
      for (size_t i=1; i<size; i++) {
        prev_key += decodeUnsigned<K>(start);
        cur_val = decodeUnsigned<V>(start);
        f(Entry::to_entry(prev_key, cur_val));
      }
    }

    template <class F>
    static inline bool decode_cond(uint8_t* bytes, uint32_t size, const F& f) {
      uint8_t* start = bytes;

      K prev_key = *((K*)start);
      start += sizeof(K);
      V cur_val = decodeUnsigned<V>(start);

      if (!f(Entry::to_entry(prev_key, cur_val))) return false;
      for (uint32_t i=1; i<size; i++) {
        prev_key += decodeUnsigned<K>(start);
        cur_val = decodeUnsigned<V>(start);
        if (!f(Entry::to_entry(prev_key, cur_val))) return false;
      }
      return true;
    }

    // F: ET -> K
    template <class F, class Comp, class K>
    static inline std::optional<ET> find(uint8_t* bytes, size_t size,
          const F& f, const Comp& comp, const K& k) {
      std::optional<ET> ret;
      auto test = [&] (const ET& et) -> bool {
        if (!(comp(k, f(et)) || comp(f(et), k))) {
          ret = et;
          return false;
        }
        return true;
      };
      decode_cond(bytes, size, test);
      return ret;
    }

    static inline void destroy(uint8_t* bytes, size_t size) { }
  };
};

}  // namespace cpam
