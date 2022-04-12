import os, glob
import sys
import signal
import time
import subprocess
import statistics
from pathlib import Path

graphs = glob.glob("../inputs/*.adj")

binaries = [("cpam", "CPAM")]

applications = ["BFS", "MIS", "BC"]

options = ["", "-flatsnap"]

binary_dir = "../algorithms/"
num_rounds = 21

Path("outputs/").mkdir(exist_ok=True)

def get_raw_output_file(name, graph, binary):
  graphname = graph.split("/")[-1].rstrip(".adj")
  return ("outputs/" + graphname + "-" + name + ".raw_output")

def get_intermediate_output_file(name, graph, binary):
  graphname = graph.split("/")[-1].rstrip(".adj")
  return ("outputs/" + graphname + "-" + name + ".int_output")

def get_output_file(name, graph, binary):
  graphname = graph.split("/")[-1].rstrip(".adj")
  return ("outputs/" + graphname + "-" + name + ".output")
