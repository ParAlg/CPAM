#!/bin/bash
# echo "BASELINE COMPARISON"
make
P=/ssd1/data/bigann
# ./neighbors -R 64 -L 128 -k 10 -Q 250 -q $P/  $P/

./neighbors -R 64 -L 128 -k 10 -Q 250 -q $P/query.public.10K.u8bin -c $P/bigann-100M -res $P/bigann_vamana.csv -t uint8 -D Euclidian $P/base.1B.u8bin.crop_nb_100000000
