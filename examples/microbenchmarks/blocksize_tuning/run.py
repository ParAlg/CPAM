import subprocess

import argparse
import settings

# Run tests based on settings and generate raw_output files.

parser = argparse.ArgumentParser(description='Run microbenchmarks')
parser.add_argument('--small', action="store_false", help='use small-size inputs')

args = parser.parse_args()
use_large = True
tests = settings.large_tests
if (args.small != None):
  use_large = False
  tests = settings.small_tests

def no_numa():
  cmd = "numactl --show"
  output = ""
  try:
    output = str(subprocess.check_output(cmd, shell=True))
  except subprocess.CalledProcessError as e:
    output = str(e.output)
  # kind of brittle...
  return "No NUMA" in output

header =  "" if no_numa() else " numactl -i all "

def run_cmd(name, binary, mode, test_num, test_name, n, m):
  raw_output_file = settings.get_raw_output_file(name, test_name, n, m)
  rounds = 1 if settings.is_sequential(binary) else settings.num_rounds
  cmd = (binary + " -s -r " + str(rounds) + " -n " + str(n) + " -m " + str(m) + " " + str(test_num) + " > " + raw_output_file)

  print(cmd)
  subprocess.call(cmd, shell=True)

for (_name, _binary) in settings.binaries:
  for block_size in settings.block_sizes:
    binary = _binary + "-" + str(block_size)
    name = _name + "-" + str(block_size)
    the_binary = settings.binary_dir + binary
    for (test, n, m, test_name) in tests:
      run_cmd(name, header + the_binary, name, test, test_name, n, m)
