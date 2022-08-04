#!/bin/bash
make clean all
rm outFileC
./neighborsCycle -R 100 -L 125 -k 200 -Q 250 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFileC /ssd1/ANN/sift/sift1M.bvecs
cd ../bench
./neighborsCheck -r [1,2,5,10,50,100,125,150,175,200] /ssd1/ANN/sift/idx_1M.ivecs ../dynamicScripts/outFileC
# cd ../dynamicScripts
# rm outFileS
#neighbors command
# cd ../bench
#check command