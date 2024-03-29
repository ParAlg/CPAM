// This code is part of the Problem Based Benchmark Suite (PBBS)
// Copyright (c) 2011 Guy Blelloch and the PBBS team
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <iostream>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include "parlay/parallel.h"
#include "parlay/primitives.h"
#include "parse_command_line.h"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

enum class ptr_mapped_src{
  MEM, DISK
};

namespace detail{

template<typename T, ptr_mapped_src Src>
class ptr_mapped_impl
{
  T *ptr_raw;
public:
  ptr_mapped_impl(){
  }

  ptr_mapped_impl(T *p) : ptr_raw(p){
  }

  template<typename U, ptr_mapped_src SrcOther>
  ptr_mapped_impl(const ptr_mapped_impl<U,SrcOther> &ptr) :
    ptr_raw(ptr.get()){
    static_assert(std::is_convertible_v<U*,T*>);
  }

  ptr_mapped_impl& operator=(T *p){
    ptr_raw = p;
    return *this;
  }

  template<typename U, ptr_mapped_src SrcOther>
  ptr_mapped_impl& operator=(const ptr_mapped_impl<U,SrcOther> &ptr){
    static_assert(std::is_convertible_v<U*,T*>);
    ptr_raw = ptr.get();
  }

  T* get() const{
    return ptr_raw;
  }

  operator T*() const{
    return get();
  }

  // For simplicity, we only keep the least methods to satisfy the requirements of LegacyIterator

  T& operator*() const{
    return *get();
  }

  ptr_mapped_impl& operator++(){
    ++ptr_raw;
    return *this;
  }
};

} // namespace detail

template<typename T, ptr_mapped_src Src>
using ptr_mapped = std::conditional_t<Src==ptr_mapped_src::MEM, T*, detail::ptr_mapped_impl<T,Src>>;

template<typename T, ptr_mapped_src Src>
struct std::iterator_traits<detail::ptr_mapped_impl<T,Src>>
{
  using difference_type = std::ptrdiff_t;
  using value_type = std::remove_cv_t<T>;
  using pointer = T*;
  using reference = T&;
  using iterator_category = void;
};

// *************************************************************
//  SOME DEFINITIONS
// *************************************************************


// *************************************************************
// Parsing code (should move to common?)
// *************************************************************

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

template<typename T, typename Conv>
auto parse_vecs(const char* filename, Conv converter)
{
  const auto [fileptr, length] = mmapStringFromFile(filename);

  // Each vector is 4 + sizeof(T)*d bytes.
  // * first 4 bytes encode the dimension (as an uint32_t)
  // * next d values are T-type variables representing vector components
  // See http://corpus-texmex.irisa.fr/ for more details.

  const uint32_t d = *((const uint32_t*)fileptr);
  std::cout << "Dimension = " << d << std::endl;

  const size_t vector_size = sizeof(d) + sizeof(T)*d;
  const size_t num_vectors = length / vector_size;
  // std::cout << "Num vectors = " << num_vectors << std::endl;

  typedef ptr_mapped<T,ptr_mapped_src::DISK> type_ptr;
  parlay::sequence<decltype(converter(0,type_ptr(nullptr),type_ptr(nullptr)))> points(num_vectors);

  parlay::parallel_for(0, num_vectors, [&] (size_t i) {
    const size_t offset_in_bytes = vector_size*i + sizeof(d);  // skip dimension
    const T* begin = (const T*)(fileptr + offset_in_bytes);
    const T* end = begin + d;
    points[i] = converter(i, type_ptr(const_cast<T*>(begin)), type_ptr(const_cast<T*>(end)));
  });

  return std::make_pair(points,d);
}

/*
auto parse_fvecs(const char* filename)
{
  return parse_vecs<float>(filename, [](size_t id, auto begin, auto end){
    typedef typename std::iterator_traits<decltype(begin)>::value_type type_elem;
    static_assert(std::is_same_v<decltype(begin),ptr_mapped<type_elem,ptr_mapped_src::DISK>>);

    Tvec_point<type_elem> point;
    point.id = id;
    point.coordinates = parlay::make_slice(begin.get(), end.get());
    return point;
  }).first;
}

auto parse_ivecs(const char* filename)
{
  return parse_vecs<float>(filename, [](size_t id, auto begin, auto end){
    typedef typename std::iterator_traits<decltype(begin)>::value_type type_elem;
    static_assert(std::is_same_v<decltype(begin),ptr_mapped<type_elem,ptr_mapped_src::DISK>>);

    ivec_point point;
    point.id = id;
    point.coordinates = parlay::make_slice(begin.get(), end.get());
    return point;
  }).first;
}
*/