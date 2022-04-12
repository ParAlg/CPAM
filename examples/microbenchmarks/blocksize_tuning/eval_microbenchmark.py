import subprocess
import statistics
import argparse

import settings

parser = argparse.ArgumentParser(description='Run microbenchmarks')
parser.add_argument('--small', action="store_false", help='use small-size inputs')

args = parser.parse_args()
use_large = True
tests = settings.large_tests
if (args.small != None):
  use_large = False
  tests = settings.small_tests


def log_results(name, test_name, n, m):
  raw_output_file = settings.get_raw_output_file(name, test_name, n, m)
  int_output_file = settings.get_intermediate_output_file(name, test_name, n, m)
  cmd = ("grep -r \"RESULT\" " + raw_output_file + " > " + int_output_file)

  print(cmd)
  subprocess.call(cmd, shell=True)

  output_file = settings.get_output_file(name, test_name, n, m)
  with open(int_output_file) as f:
    lines = f.readlines()
    times = []
    for line in lines:
      l = float(line.split("=")[2].split("\t")[0])
      times.append(l)
    median = round(statistics.median(times), 5)

    output_str = name + "," + test_name + "," + str(median) + "\n"
    print(output_file)
    with open(output_file, "w") as of:
      of.write(output_str)


for (_name, _binary) in settings.binaries:
  for block_size in settings.block_sizes:
    binary = _binary + "-" + str(block_size)
    name = _name + "-" + str(block_size)
    the_binary = settings.binary_dir + binary
    for (test, n, m, test_name) in tests:
      log_results(name, test_name, n, m)

cmd = ("cat outputs/*.output  > benchmark_results.txt")
print(cmd)
subprocess.call(cmd, shell=True)
