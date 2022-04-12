package(
  default_visibility = ["//visibility:public"],
)

cc_library(
  name = "augmented_map",
  hdrs = ["augmented_map.h"],
  deps = []
)

cc_library(
  name = "augmented_node",
  hdrs = ["augmented_node.h"],
  deps = [
  ":basic_node",
  ]
)

cc_library(
  name = "augmented_ops",
  hdrs = ["augmented_ops.h"],
  deps = [
  ":map_ops",
  ":utils",
  ]
)

cc_library(
  name = "balance_utils",
  hdrs = ["balance_utils.h"],
  deps = []
)

cc_library(
  name = "basic_node",
  hdrs = ["basic_node.h"],
  deps = [
  ":basic_node_helpers",
  ":byte_encode",
  ":compression",
  ":utils",
  "//parlaylib/include/parlay:alloc",
  ]
)

cc_library(
  name = "basic_node_helpers",
  hdrs = ["basic_node_helpers.h"],
  deps = [
  ":utils",
  ]
)

cc_library(
  name = "build",
  hdrs = ["build.h"],
  deps = [
  ":get_time",
  "//parlaylib/include/parlay:primitives",
  ]
)

cc_library(
  name = "byte_encode",
  hdrs = ["byte_encode.h"],
  deps = []
)

cc_library(
  name = "compression",
  hdrs = ["compression.h"],
  deps = []
)

cc_library(
  name = "gc",
  hdrs = ["gc.h"],
  deps = [
  ":basic_node",
  ":utils",
  ]
)

cc_library(
  name = "get_time",
  hdrs = ["get_time.h"],
  deps = []
)

cc_library(
  name = "map",
  hdrs = ["map.h"],
  deps = []
)

cc_library(
  name = "map_ops",
  hdrs = ["map_ops.h"],
  deps = [
  ":utils",
  "//parlaylib/include/parlay:primitives",
  "//parlaylib/include/parlay/internal:binary_search",
  ]
)

cc_library(
  name = "cpam",
  hdrs = ["cpam.h"],
  deps = [
  ":basic_node",
  ":augmented_node",
  ":utils",
  ":gc",
  ":weight_balanced_tree",
  ":sequence_ops",
  ":map_ops",
  ":augmented_ops",
  ":build",
  ":map",
  ":augmented_map",
  "//parlaylib/include/parlay:utilities",
  ]
)

cc_library(
  name = "parse_command_line",
  hdrs = ["parse_command_line.h"],
  deps = []
)

cc_library(
  name = "sequence_ops",
  hdrs = ["sequence_ops.h"],
  deps = [
  ":gc",
  ":utils",
  ]
)

cc_library(
  name = "utils",
  hdrs = ["utils.h"],
  deps = []
)

cc_library(
  name = "weight_balanced_tree",
  hdrs = ["weight_balanced_tree.h"],
  deps = [
  ":balance_utils",
  ]
)
