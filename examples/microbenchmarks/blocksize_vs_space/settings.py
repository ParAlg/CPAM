import os
import sys
import signal
import time
import statistics
from pathlib import Path

binaries = [
  ("cpam-na", "testParallel-CPAM-NA"),
  ("cpam-na-diff", "testParallel-CPAM-NA-Diff"),
  ("cpam", "testParallel-CPAM"),
  ("cpam-diff", "testParallel-CPAM-Diff"),
]

pam_binaries = [
  ("pam", "testParallel-PAM"),
  ("pam-na", "testParallel-PAM-NA"),
]

block_sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]

large_tests = [
  (35, 1e8, 1e8,  "size_in_bytes"),   # == Large Tests == use for large runs.
]

small_tests = [
  (35, 1e5, 1e5,  "size_in_bytes"),  # == Small Tests == use for small runs.
]

binary_dir = "../"
output_dir = "outputs/"
num_rounds = 0   # warmup run is sufficient; size is deterministic

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
