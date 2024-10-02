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
  name = "avl_tree",
  hdrs = ["avl_tree.h"],
  deps = [
  ":balance_utils",
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
  "@parlaylib//parlay:alloc",
  ]
)

cc_library(
  name = "build",
  hdrs = ["build.h"],
  deps = [
  ":get_time",
  "@parlaylib//parlay:primitives",
  ]
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
  "@parlaylib//parlay:primitives",
  "@parlaylib//parlay/internal:binary_search",
  ]
)

cc_library(
  name = "pam",
  hdrs = ["pam.h"],
  deps = [
  ":basic_node",
  ":augmented_node",
  ":utils",
  ":gc",
  ":avl_tree",
  ":red_black_tree",
  ":weight_balanced_tree",
  ":treap",
  ":sequence_ops",
  ":map_ops",
  ":augmented_ops",
  ":build",
  ":map",
  ":augmented_map",
  "@parlaylib//parlay:utilities",
  ]
)

cc_library(
  name = "parse_command_line",
  hdrs = ["parse_command_line.h"],
  deps = []
)

cc_library(
  name = "red_black_tree",
  hdrs = ["red_black_tree.h"],
  deps = [
  ":balance_utils",
  ]
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
  name = "treap",
  hdrs = ["treap.h"],
  deps = [
  ":balance_utils",
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
