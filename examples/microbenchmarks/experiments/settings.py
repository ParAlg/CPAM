import os
import sys
import signal
import time
import statistics
from pathlib import Path

unaug_binaries = [
  ("cpam-na", "testParallel-CPAM-NA"),
  ("cpam-na-diff", "testParallel-CPAM-NA-Diff"),
  ("pam-na", "testParallel-PAM-NA"),
  ("cpam-na-seq", "testParallel-CPAM-NA-Seq"),
  ("cpam-na-diff-seq", "testParallel-CPAM-NA-Diff-Seq"),
  ("pam-na-seq", "testParallel-PAM-NA-Seq"),
]

aug_binaries = [
  ("cpam", "testParallel-CPAM"),
  ("cpam-diff", "testParallel-CPAM-Diff"),
  ("pam", "testParallel-PAM"),
  ("cpam-seq", "testParallel-CPAM-Seq"),
  ("cpam-diff-seq", "testParallel-CPAM-Diff-Seq"),
  ("pam-seq", "testParallel-PAM-Seq"),
]

# LARGE tests
large_unaug_tests = [
  (0,  1e8, 1e5,  "union"),
  (0,  1e8, 1e8,  "union"),
  (4,  1e8, 1e8,  "build"),
  (1,  1e8, 1e8,  "intersect"),
  (9,  1e8, 1e8,  "difference"),
  (33, 1e8, 1e8,  "map"),
  (34, 1e8, 1e8,  "reduce"),
  (5,  1e8, 1e8,  "filter"),
  (14, 1e8, 1e8, "multi_insert"),
  (12, 1e8, 1e8, "find"),
  (2,  1e8, 1e6,  "insert"),
  (18, 1e8, 1e6, "range"),
]

large_aug_tests = [
  (0, 1e8, 1e8,  "union"),
  (4, 1e8, 1e8,  "build"),
  (18, 1e8, 1e6, "aug_range"),
  (18, 1e8, 1e6, "aug_filter"),
]

# SMALL tests. Uncomment to test out the code / run the tests faster.
# Note that the results will be different on the smaller tests.
small_unaug_tests = [
  (0,  1e6, 1e4,  "union"),
  (0,  1e6, 1e6,  "union"),
  (4,  1e6, 1e6,  "build"),
  (1,  1e6, 1e6,  "intersect"),
  (9,  1e6, 1e6,  "difference"),
  (33,  1e6, 1e6,  "map"),
  (34,  1e6, 1e6,  "reduce"),
  (5,  1e6, 1e6,  "filter"),
  (14, 1e6, 1e6, "multi_insert"),
  (12, 1e6, 1e4, "find"),
  (2,  1e6, 1e4,  "insert"),
  (18, 1e6, 1e4, "range"),
]

small_aug_tests = [
  (4,  1e6, 1e6,  "build"),
  (0,  1e6, 1e6,  "union"),
  (18, 1e6, 1e4, "aug_range"),
  (18, 1e6, 1e4, "aug_filter"),
]

binary_dir = "../"
output_dir = "outputs/"
num_rounds = 7

Path(output_dir).mkdir(exist_ok=True)

def is_sequential(binary_name):
  return "Seq" in binary_name

def get_raw_output_file(name, test_name, n, m):
  the_name = name + "-" + test_name + "-" + str(int(n)) + "-" + str(int(m))
  return ("outputs/" + the_name + ".raw_output")

def get_intermediate_output_file(name, test_name, n, m):
  the_name = name + "-" + test_name + "-" + str(int(n)) + "-" + str(int(m))
  return ("outputs/" + the_name + ".int_output")

def get_output_file(name, test_name, n, m):
  the_name = name + "-" + test_name + "-" + str(int(n)) + "-" + str(int(m))
  return ("outputs/" + the_name + ".output")
