#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

export THREADS=`nproc --all`

if [[ -z "${SMALL}" ]]; then
small_flag=
else
small_flag="--small"
fi

printf "${GREEN}Running Microbenchmarks${NC}\n"
printf "(Results in Figure 10)\n"
cd ../
make -j $THREADS
cd blocksize_vs_space

python3 run.py $small_flag
python3 eval_microbenchmark.py $small_flag

echo
printf "${GREEN}RESULTS${NC}\n"
cat benchmark_results.txt