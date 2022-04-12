import os
import sys
import signal
import time
import subprocess
import statistics
from pathlib import Path

import settings

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
  cmd = (binary + " -src 10012 -s -m -rounds " + str(settings.num_rounds) + " " + graph + " > " + raw_output_file)

  subprocess.call(cmd, shell=True)

for graph in settings.graphs:
  for application in settings.applications:
    for (name, binary) in settings.binaries:
      for option in settings.options:
        the_binary = settings.binary_dir + application + "-" + binary
        if (option != ""):
          the_binary = the_binary + " " + option
        the_name = application + "-" + name + option
        run_cmd(the_name, graph, header + the_binary, application, name)

for graph in settings.graphs:
  for (name, binary) in settings.binaries:
    for option in settings.options:
      the_binary = settings.binary_dir + "Flatsnap" + "-" + binary
      if (option != ""):
        the_binary = the_binary + " " + option
      the_name = "FLATSNAP" + "-" + name + option
      run_cmd(the_name, graph, header + the_binary, "FLATSNAP", name)
