licenses(["notice"])

cc_library(
    name = "alloc",
    hdrs = ["alloc.h"],
    deps = [
        ":parallel",
        ":utilities",
        "//parlay/internal:block_allocator",
        "//parlay/internal:concurrent_stack",
        "//parlay/internal:memory_size",
    ],
)

cc_library(
    name = "delayed_sequence",
    hdrs = ["delayed_sequence.h"],
)

cc_library(
    name = "hash_table",
    hdrs = ["hash_table.h"],
    deps = [
        ":delayed_sequence",
        ":monoid",
        ":parallel",
        ":sequence",
        ":slice",
        ":utilities",
        "//parlay/internal:sequence_ops",
    ],
)

cc_library(
    name = "io",
    hdrs = ["io.h"],
    deps = [
        ":primitives",
        ":sequence",
        ":slice",
        "//parlay/internal:file_map",
    ],
)

cc_library(
    name = "monoid",
    hdrs = ["monoid.h"],
)

cc_library(
    name = "parallel",
    hdrs = ["parallel.h"],
    deps = [
        ":scheduler",
        "//parlay/internal/scheduler_plugins:cilk",
        "//parlay/internal/scheduler_plugins:omp",
        "//parlay/internal/scheduler_plugins:sequential",
    ],
)

cc_library(
    name = "primitives",
    hdrs = ["primitives.h"],
    deps = [
        ":delayed_sequence",
        ":monoid",
        ":parallel",
        ":range",
        ":sequence",
        ":slice",
        ":utilities",
        "//parlay/internal:block_delayed",
        "//parlay/internal:collect_reduce",
        "//parlay/internal:integer_sort",
        "//parlay/internal:merge",
        "//parlay/internal:merge_sort",
        "//parlay/internal:sample_sort",
        "//parlay/internal:sequence_ops",
    ],
)

cc_library(
    name = "random",
    hdrs = ["random.h"],
    deps = [
        ":delayed_sequence",
        ":parallel",
        ":range",
        ":sequence",
        ":slice",
        ":utilities",
        "//parlay/internal:counting_sort",
    ],
)

cc_library(
    name = "range",
    hdrs = ["range.h"],
)

cc_library(
    name = "scheduler",
    hdrs = ["scheduler.h"],
    deps = [
        "//parlay/internal:work_stealing_job",
    ],
)

cc_library(
    name = "sequence",
    hdrs = ["sequence.h"],
    deps = [
        ":alloc",
        ":parallel",
        ":range",
        ":slice",
        ":type_traits",
        ":utilities",
        "//parlay/internal:debug_uninitialized",
        "//parlay/internal:sequence_base",
    ],
)

cc_library(
    name = "slice",
    hdrs = ["slice.h"],
    deps = [
        ":range",
    ],
)

cc_library(
    name = "type_traits",
    hdrs = ["type_traits.h"],
)

cc_library(
    name = "utilities",
    hdrs = ["utilities.h"],
    deps = [
        ":parallel",
        ":type_traits",
        "//parlay/internal:debug_uninitialized",
    ],
)

package(
    default_visibility = ["//visibility:public"],
)
