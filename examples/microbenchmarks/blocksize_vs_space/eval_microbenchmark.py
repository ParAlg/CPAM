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

output_file = "benchmark_results.txt"

cmd = ("rm benchmark_results.txt")
subprocess.call(cmd, shell=True)

def log_results(name, test_name, n, m, of):
  raw_output_file = settings.get_raw_output_file(name, test_name, n, m)
  int_output_file = settings.get_intermediate_output_file(name, test_name, n, m)
  cmd = ("grep -r \"size in bytes\" " + raw_output_file + " > " + int_output_file)

  subprocess.call(cmd, shell=True)

  # output_file = settings.get_output_file(name, test_name, n, m)
  # print(int_output_file)
  with open(int_output_file) as f:
    lines = f.readlines()
    times = []
    for line in lines:
      l = float(line.split("=")[-1])
      times.append(l)
    median = round(statistics.median(times), 5)

    output_str = name + "," + test_name + "," + str(median) + "\n"
    of.write(output_str)

with open(output_file, "w") as of:
  for (_name, _binary) in settings.binaries:
    for block_size in settings.block_sizes:
      binary = _binary + "-" + str(block_size)
      name = _name + "-" + str(block_size)
      the_binary = settings.binary_dir + binary
      for (test, n, m, test_name) in tests:
        log_results(name, test_name, n, m, of)

  for (_name, _binary) in settings.pam_binaries:
    binary = _binary
    name = _name
    the_binary = settings.binary_dir + binary
    for (test, n, m, test_name) in tests:
      log_results(name, test_name, n, m, of)


print("Wrote results to " + output_file)

#cmd = ("cat outputs/*.output  > benchmark_results.txt")
#print(cmd)
#subprocess.call(cmd, shell=True)

