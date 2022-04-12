import os
import sys
import signal
import time
import subprocess
import statistics
from pathlib import Path

import settings

data = {}

def log_results(name, graph, binary, app_name, mode):
  raw_output_file = settings.get_raw_output_file(name, graph, binary)
  int_output_file = settings.get_intermediate_output_file(name, graph, binary)
  cmd = ("grep -r \"### Running Time: \" " + raw_output_file + " | sed 's/### Running Time: //' > " + int_output_file)
  if ("aspen" in binary):
    cmd = ("grep -r \"RESULT\ttest=\" " + raw_output_file + " | sed 's/RESULT\t//' > " + int_output_file)
  graphname = graph.split("/")[-1].rstrip(".adj")

  subprocess.call(cmd, shell=True)

  output_file = settings.get_output_file(name, graph, binary)
  with open(int_output_file) as f:
    lines = f.readlines()
    times = []
    for line in lines:
      l = ""
      if ("aspen" in binary):
        l = float(line.split("\t")[1].split("=")[1])
      else:
        l = float(line.rstrip())
      times.append(l)
    median = statistics.median(times)

    if graph not in data:
      data[graph] = {}
    if (app_name == "FLATSNAP"):
      data[graph][mode + "-FS"] = median

    using_fs = "-flatsnap" in name
    if (using_fs):
      data[graph][mode + "-" + app_name + "-FS"] = median
    else:
      data[graph][mode + "-" + app_name] = median

    output_str = app_name + "," + graphname + "," + mode + "," + str(median) + "\n"
    print(output_str)
#    print(output_file)
#    with open(output_file, "w") as of:
#      of.write(output_str)

all_apps = settings.applications + ["FLATSNAP"]

for graph in settings.graphs:
  for application in all_apps:
    for (name, binary) in settings.binaries:
      for option in settings.options:
        the_binary = settings.binary_dir + application + "-" + binary
        the_name = application + "-" + name + option
        log_results(the_name, graph, " numactl -i all " + the_binary, application, name)

aspen_binaries = [("aspen", "run_static_algorithm")]
aspen_binary_dir = "../other_systems/aspen/"

for graph in settings.graphs:
  for application in all_apps:
    for (name, binary) in aspen_binaries:
      for option in settings.options:
        the_binary = aspen_binary_dir + binary + " -t " + application
        if (option != ""):
          the_binary = the_binary + " " + option
        the_name = application + "-" + name + option
        log_results(the_name, graph, " numactl -i all " + the_binary, application, name)

output_file = "benchmark_results.txt"

def mystr(val):
  return format(val, ".4f")

with open(output_file, "w") as f:
  for app in settings.applications:
    for graph in settings.graphs:
      cpam = "cpam"
      aspen_fs = data[graph]["aspen-" + app + "-FS"]
      aspen_fs_time = data[graph]["aspen" + "-FS"]

      no_fs = data[graph][cpam + "-" + app]
      with_fs = data[graph][cpam + "-" + app + "-FS"]
      no_fs_over_fs = no_fs / with_fs

      our_fs_time = data[graph][cpam + "-FS"]

      aspen_over_ours = aspen_fs / with_fs

      output_str = "Application = " + app + ", Graph = " + graph + ", Aspen FS = " + mystr(aspen_fs) + ", Aspen FS time = " + mystr(aspen_fs_time) + ", Ours No-FS = " + mystr(no_fs) + ", Ours FS = " + mystr(with_fs) + ", No-FS/FS = " + mystr(no_fs_over_fs) + ", Ours FS time = " + mystr(our_fs_time) + ", Aspen/Ours = " + mystr(aspen_over_ours) + "\n"
      f.write(output_str)
