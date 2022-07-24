import subprocess
import statistics

import argparse
import settings

results_file = "microbenchmark_results.txt"

tab = {}

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

def log_results(name, test_name, n, m, seq):
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

    trunc_name = name.rstrip("-seq")
    if (trunc_name not in tab):
      tab[trunc_name] = {}
    if (test_name not in tab[trunc_name]):
      tab[trunc_name][test_name] = {}
    mode_str = "Seq" if seq else "Par"
    tab[trunc_name][test_name][mode_str] = median

    output_str = name + "," + test_name + "," + str(median) + "\n"
    print(output_file)
    with open(output_file, "w") as of:
      of.write(output_str)

# To emit the .output files:
for (name, binary) in settings.unaug_binaries:
  for (test, n, m, test_name) in unaug_tests:
    if (settings.is_sequential(binary)):
      log_results(name, test_name, n, m, True)
    else:
      log_results(name, test_name, n, m, False)

for (name, binary) in settings.aug_binaries:
  for (test, n, m, test_name) in aug_tests:
    if (settings.is_sequential(binary)):
      log_results(name, test_name, n, m, True)
    else:
      log_results(name, test_name, n, m, False)

# Compute speedup numbers and print microbench_results.

strs = []
for (name, binary) in settings.unaug_binaries:
  for (test, n, m, test_name) in unaug_tests:
    if (settings.is_sequential(binary)):  # symmetry break
      trunc_name = name.rstrip("-seq")
      seq = tab[trunc_name][test_name]["Seq"]
      par = tab[trunc_name][test_name]["Par"]
      spd = format(seq / par, ".2f")
      output_str = "system = " + trunc_name.replace("-na-", "-").replace("-na", "") + ", mode = NoAug, name = " + test_name + ", n = " + str(n) + ", m = " + str(m) + ", sequential_time = " + str(seq) + ", parallel time = " + str(par) + ", speedup = " + str(spd) + "\n"
      strs.append(output_str)

for (name, binary) in settings.aug_binaries:
  for (test, n, m, test_name) in aug_tests:
    if (settings.is_sequential(binary)):  # symmetry break
      trunc_name = name.rstrip("-seq")
      seq = tab[trunc_name][test_name]["Seq"]
      par = tab[trunc_name][test_name]["Par"]
      spd = format(seq / par, ".2f")
      output_str = "system = " + trunc_name.replace("-na-", "-").replace("-na", "") + ", mode = Aug, name = " + test_name + ", n = " + str(n) + ", m = " + str(m) + ", sequential_time = " + str(seq) + ", parallel time = " + str(par) + ", speedup = " + str(spd) + "\n"
      strs.append(output_str)

strs.sort()
with open(results_file, "w") as f:
  for l in strs:
    f.write(l)


