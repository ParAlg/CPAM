#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

export THREADS=`nproc --all`

if [[ -z "${SMALL}" ]]; then
n=100000000
q=100000000
else
n=1000000
q=1000000
fi

r=5

header="numactl -i all"
numa_exists="$(numactl --show 2>&1)"
if [[ $numa_exists == *"No NUMA"* ]]; then
    header=""
fi

printf "${GREEN}test interval tree (CPAM)${NC}\n"
printf "(Results in Table 3, row 4)\n"
make -j $THREADS

printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./intervalTree -r $r -n $n -q $q 7 | grep "RESULT" | tee res.txt

echo
printf "${BLUE}sequential run${NC}\n"
export PARLAY_NUM_THREADS=1
./intervalTree_seq -r $r -n $n -q $q 3 | grep "RESULT" | tee -a res.txt

echo
printf "${GREEN}test interval tree (PAM)${NC}\n"
printf "(Results in Table 3, row 5)\n"
make -j $THREADS

printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./intervalTree_pam -r $r -n $n -q $q 7 | grep "RESULT" | tee res_pam.txt

echo
printf "${BLUE}sequential run${NC}\n"
export PARLAY_NUM_THREADS=1
./intervalTree_pam_seq -r $r -n $n -q $q 3 | grep "RESULT" | tee -a res_pam.txt

python3 comp.py