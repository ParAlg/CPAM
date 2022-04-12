import os, glob

files = glob.glob("./*.output")
print(files)

output = "benchmark_output.txt"

def gib(val):
  return float(val) / (1024**3)

def my_str(val):
  return format(val, ".4f")

with open(output, "w") as out:
  for fname in files:
    with open(fname) as f:
      l1 = f.readline()
      l2 = f.readline()
      l3 = f.readline()
      l4 = f.readline()
      l5 = f.readline()
      l6 = f.readline()

      pam_pam = gib(float(l1.split(",")[4]))
      pam_cpam = gib(float(l2.split(",")[4]))
      cpam_cpam = gib(float(l3.split(",")[4]))
      cpam_cpam_diff = gib(float(l4.split(",")[4]))
      aspen = float(l5.split("=")[1].rstrip())
      gbbs = gib(float(l6.split("\t")[0]))

      output_str = "Graph = " + fname.lstrip("./").rstrip(".output") + ", GBBS = (" + my_str(gbbs) + ", " + my_str(gbbs/gbbs) + "), PaC-tree (Diff) = (" + my_str(cpam_cpam_diff) + ", " + my_str(cpam_cpam_diff/gbbs) + "), PaC-tree = (" + my_str(cpam_cpam) + ", " + my_str(cpam_cpam/gbbs) + "), PAM = (" + my_str(pam_pam) + ", " + my_str(pam_pam/gbbs) + "), Aspen = (" + my_str(aspen) + ", " + my_str(aspen/gbbs) + ")\n"
      out.write(output_str)
