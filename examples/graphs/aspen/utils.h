#pragma once

#include "parlay/primitives.h"

namespace aspen {

  struct empty { };

  using vertex_id = uint32_t;
  using edge_id = size_t;

  // Writes the list of indices `i` where `Fl[i] == true` to range `Out`.
  template <class Bool_Seq, class Out_Seq>
  size_t pack_index_out(Bool_Seq const &Fl, Out_Seq&& Out) {
    using Idx_Type = typename std::remove_reference<Out_Seq>::type::value_type;
    auto identity = [] (size_t i) {return (Idx_Type) i;};
    return parlay::pack_into(
	parlay::delayed_seq<Idx_Type>(Fl.size(),identity),
	Fl,
	std::forward<Out_Seq>(Out));
  }

  template <class Idx_Type, class D, class F>
  inline parlay::sequence<std::tuple<Idx_Type, D> > pack_index_and_data(
      F& f, size_t size) {
    auto id_seq = parlay::delayed_seq<std::tuple<Idx_Type, D> >(size,  [&](size_t i) {
      return std::make_tuple((Idx_Type)i, std::get<1>(f[i]));
    });
    auto flgs_seq = parlay::delayed_seq<bool>(size, [&](size_t i) { return std::get<0>(f[i]); });

    return parlay::pack(id_seq, flgs_seq);
  }

}  // namespace aspen
