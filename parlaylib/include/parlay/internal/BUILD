cc_library(
  name = "binary_search",
  hdrs = ["binary_search.h"],
)

cc_library(
  name = "block_allocator",
  hdrs = ["block_allocator.h"],
  deps = [
    ":concurrent_stack",
    ":memory_size",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "bucket_sort",
  hdrs = ["bucket_sort.h"],
  deps = [
    "//parlay:utilities",
  ],
)

cc_library(
  name = "collect_reduce",
  hdrs = ["collect_reduce.h"],
  deps = [
    "//parlay:utilities",
  ],
)

cc_library(
  name = "concurrent_stack",
  hdrs = ["concurrent_stack.h"],
  deps = [
    "//parlay:utilities",
  ],
)

cc_library(
  name = "counting_sort",
  hdrs = ["counting_sort.h"],
  deps = [
    ":sequence_ops",
    ":transpose",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "debug_uninitialized",
  hdrs = ["debug_uninitialized.h"],
)

cc_library(
  name = "file_map",
  hdrs = ["file_map.h"],
  deps = [
    "//parlay/internal/posix:file_map_impl_posix",
    "//parlay/internal/windows:file_map_impl_windows",
  ],
)

cc_library(
  name = "integer_sort",
  hdrs = ["integer_sort.h"],
  deps = [
    ":counting_sort",
    ":sequence_ops",
    ":quicksort",
    ":uninitialized_sequence",
    "//parlay:delayed_sequence",
    "//parlay:slice",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "memory_size",
  hdrs = ["memory_size.h"],
)

cc_library(
  name = "merge",
  hdrs = ["merge.h"],
  deps = [
    "//parlay:sequence",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "merge_sort",
  hdrs = ["merge_sort.h"],
  deps = [
    ":merge",
    ":quicksort",
    "//parlay:sequence",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "quicksort",
  hdrs = ["quicksort.h"],
  deps = [
    ":uninitialized_storage",
    ":uninitialized_sequence",
    ":sequence_ops",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "sample_sort",
  hdrs = ["sample_sort.h"],
  deps = [
    ":bucket_sort",
    ":quicksort",
    ":sequence_ops",
    ":transpose",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "sequence_base",
  hdrs = ["sequence_base.h"],
  deps = [
    "//parlay:parallel",
    "//parlay:type_traits",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "sequence_ops",
  hdrs = ["sequence_ops.h"],
  deps = [
    "//parlay:delayed_sequence",
    "//parlay:monoid",
    "//parlay:range",
    "//parlay:sequence",
    "//parlay:utilities",
  ],
)

cc_library(
  name = "transpose",
  hdrs = ["transpose.h"],
  deps = [
    "//parlay:utilities",
  ],
)

cc_library(
  name = "uninitialized_sequence",
  hdrs = ["uninitialized_sequence.h"],
  deps = [
    ":debug_uninitialized",
  ],
)

cc_library(
  name = "uninitialized_storage",
  hdrs = ["uninitialized_storage.h"],
  deps = [
    ":debug_uninitialized",
  ],
)

cc_library(
  name = "work_stealing_job",
  hdrs = ["work_stealing_job.h"],
)

package(
  default_visibility = ["//visibility:public"],
)
