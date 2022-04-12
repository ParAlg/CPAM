#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

export THREADS=`nproc --all`

printf "${GREEN}Running Microbenchmarks${NC}\n"
printf "(Results in Table 2 (Microbenchmark results))\n"
make -j $THREADS

if [[ -z "${SMALL}" ]]; then
small_flag=
else
small_flag="--small"
fi

cd experiments
python3 run_microbenchmark.py $small_flag
python3 eval_microbenchmark.py $small_flag
cd ../

echo
printf "${GREEN}RESULTS${NC}\n"
cat experiments/microbenchmark_results.txt