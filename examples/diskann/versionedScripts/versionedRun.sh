#!/bin/bash
make clean all
./insertsAndQueries -R 64 -L 128 -k 200 -Q 400 /ssd1/ANN/sift/sift1M.bvecs