import os, glob
import sys
import signal
import time
import subprocess
import statistics
from pathlib import Path

import settings

aspen_binaries = [("aspen", "run_static_algorithm")]
aspen_binary_dir = "../other_systems/aspen/"

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

def run_cmd(name, graph, binary, app_name, mode):
  raw_output_file = settings.get_raw_output_file(name, graph, binary)
  cmd = (binary + " -src 10012 -s -m -rounds " + str(settings.num_rounds) + " -f " + graph + " > " + raw_output_file)

  subprocess.call(cmd, shell=True)

for graph in settings.graphs:
  for application in settings.applications:
    for (name, binary) in aspen_binaries:
      for option in settings.options:
        the_binary = aspen_binary_dir + binary + " -t " + application
        if (option != ""):
          the_binary = the_binary + " " + option
        the_name = application + "-" + name + option
        run_cmd(the_name, graph, header + the_binary, application, name)

for graph in settings.graphs:
  for (name, binary) in aspen_binaries:
    for option in settings.options:
      the_binary = aspen_binary_dir + binary + " -t " + "FLATSNAP"
      if (option != ""):
        the_binary = the_binary + " " + option
      the_name = "FLATSNAP" + "-" + name + option
      run_cmd(the_name, graph, header + the_binary, "FLATSNAP", name)
