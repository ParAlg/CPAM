#!/bin/bash
# echo "BASELINE COMPARISON"
cd ..
make
rm outFile
./neighbors -R 100 -L 125 -k 200 -Q 250 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFile /ssd1/ANN/sift/sift1M.bvecs
cd bench
./neighborsCheck -r [1,2,5,10,50,100,150,200] /ssd1/ANN/sift/idx_1M.ivecs ../outFile
# echo "CYCLE EXPERIMENT"
# cd ../dynamicScripts
# make clean all
# rm outFileC
# ./neighborsCycle -R 100 -L 125 -k 200 -Q 250 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFileC /ssd1/ANN/sift/sift1M.bvecs
# cd ../bench
# ./neighborsCheck -r [1,2,5,10,50,100,150,200] /ssd1/ANN/sift/idx_1M.ivecs ../dynamicScripts/outFileC
# echo "SLIDING WINDOW EXPERIMENT"
# cd ../dynamicScripts
# rm outFileS
# ./neighborsSlidingWindow -R 100 -L 125 -k 200 -Q 250 -q /ssd1/ANN/sift/bigann_query.bvecs -o outFileS /ssd1/ANN/sift/sift10M.bvecs
# cd ../bench
# ./neighborsCheck -r [1,2,5,10,50,100,150,200] /ssd1/ANN/sift/idx_1M.ivecs ../dynamicScripts/outFileS
