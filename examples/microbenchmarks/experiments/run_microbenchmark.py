import subprocess

import argparse
import settings

# Run tests based on settings and generate raw_output files.

parser = argparse.ArgumentParser(description='Run microbenchmarks')
parser.add_argument('--small', action="store_true", help='use small-size inputs')

args = parser.parse_args()
use_large = True
unaug_tests = settings.large_unaug_tests
aug_tests = settings.large_aug_tests
if args.small:
  use_large = False
  unaug_tests = settings.small_unaug_tests
  aug_tests = settings.small_aug_tests

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

for (name, binary) in settings.unaug_binaries:
  the_binary = settings.binary_dir + binary
  for (test, n, m, test_name) in unaug_tests:
    if (settings.is_sequential(binary)):
      run_cmd(name, the_binary, name, test, test_name, n, m)
    else:
      run_cmd(name, header + the_binary, name, test, test_name, n, m)

for (name, binary) in settings.aug_binaries:
  the_binary = settings.binary_dir + binary
  for (test, n, m, test_name) in aug_tests:
    if (settings.is_sequential(binary)):
      run_cmd(name, the_binary, name, test, test_name, n, m)
    else:
      run_cmd(name, header + the_binary, name, test, test_name, n, m)
