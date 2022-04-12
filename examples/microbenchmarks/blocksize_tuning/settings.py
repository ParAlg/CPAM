import os
import sys
import signal
import time
import statistics
from pathlib import Path

binaries = [
  ("cpam-na", "testParallel-CPAM-NA"),
]

block_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]

# == Large Tests== used in the paper.
large_tests = [
  (0, 1e8, 1e8,  "union"),
  (0, 1e8, 1e5,  "union-imbal"),
  (1, 1e8, 1e8,  "intersect"),
  (9, 1e8, 1e8,  "difference"),
  (4, 1e8, 1e8,  "build"),
  (5, 1e8, 1e8,  "filter"),
  (14, 1e8, 1e8, "multi_insert"),
  (12, 1e8, 1e8, "find"),
  (2, 1e8, 2e5,  "insert"),
  (18, 1e8, 1e6, "range"),
]

# == Small Tests== use these to test the scripts.
small_tests = [
  (0,  1e6, 1e6, "union"),
  (0,  1e6, 1e3, "union-imbal"),
  (1,  1e6, 1e6, "intersect"),
  (9,  1e6, 1e6, "difference"),
  (4,  1e6, 1e6, "build"),
  (5,  1e6, 1e6, "filter"),
  (14, 1e6, 1e6, "multi_insert"),
  (12, 1e6, 1e6, "find"),
  (2,  1e6, 2e4, "insert"),
  (18, 1e6, 1e4, "range"),
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
