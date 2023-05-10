#!/bin/bash
# echo "BASELINE COMPARISON"
make 
P=/ssd1/data/bigann
./neighbors -R 64 -L 128 -t uint8 -D Euclidian $P/base.1B.u8bin.crop_nb_1000000
