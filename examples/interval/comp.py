files = ["res.txt", "res_pam.txt"]
data = {
  "res.txt" : "data_cpam.txt",
  "res_pam.txt" : "data_pam.txt",
}

for fname in files:
  f = open(fname)
  line1 = f.readline()
  line2 = f.readline()
  line3 = f.readline()
  line4 = f.readline()

  items1 = line1.split(', ')
  items2 = line2.split(', ')
  items3 = line3.split(', ')
  items4 = line4.split(', ')

  build_par = float(items1[7][7:])
  build_seq = float(items3[7][7:])
  build_spd = build_seq/build_par

  query_par = float(items2[7][7:])
  query_seq = float(items4[7][7:])
  query_spd = query_seq/query_par

  size_in_gib = float(items1[8].split("=")[-1])

  fout = open(data[fname],'w')
  fout.write(items1[2]+', '+items1[3]+', ')
  fout.write("Space (GiB) = " + str(size_in_gib) + ", sequential time = " + str(build_seq) + ", parallel time = " + str(build_par) + ", speedup = " + str(build_spd)+'\n')
  fout.write(items2[2]+', '+items2[3]+', '+ items2[6] + ", ")
  fout.write("Space (GiB) = " + str(size_in_gib) + ", sequential time = " + str(query_seq) + ", parallel time = " + str(query_par) + ", speedup = " + str(query_spd)+'\n')
  fout.close()
