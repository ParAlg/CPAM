#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

export THREADS=`nproc --all`

if [[ -z "${SMALL}" ]]; then
n=100000000
q_sum=1000000
q_all=1000
else
n=100000
q_sum=100000
q_all=1000
fi

header="numactl -i all"
numa_exists="$(numactl --show 2>&1)"
if [[ $numa_exists == *"No NUMA"* ]]; then
    header=""
fi

printf "${GREEN}test range tree (CPAM)${NC}\n"
printf "(Results in Table 3, row 6)\n"
make -j $THREADS

# Run Range Sum
printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./range_test -n $n -r 5 -q $q_sum -t 1 -w 200000000 -d 1 | grep "RESULT" | tee res.txt

echo
printf "${BLUE}sequential run${NC}\n"
export PARLAY_NUM_THREADS=1
$header ./range_test_seq -n $n -r 5 -q $q_sum -t 1 -w 200000000 -d 1 | grep "RESULT" | tee -a res.txt

# Run Range All
printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./range_test -n $n -r 5 -q $q_all -t 0 -w 200000000 -d 1 | grep "RESULT" | tee -a res.txt

echo
printf "${BLUE}sequential run${NC}\n"
export PARLAY_NUM_THREADS=1
$header ./range_test_seq -n $n -r 5 -q $q_all -t 0 -w 200000000 -d 1 | grep "RESULT" | tee -a res.txt


echo
printf "${GREEN}test range tree (PAM)${NC}\n"
printf "(Results in Table 3, row 7)\n"
make -j $THREADS

# Run Range Sum
printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./range_test_pam -n $n -r 5 -q $q_sum -t 1 -w 200000000 -d 1 | grep "RESULT" | tee res_pam.txt

echo
printf "${BLUE}sequential run${NC}\n"
export PARLAY_NUM_THREADS=1
$header ./range_test_pam_seq -n $n -r 5 -q $q_sum -t 1 -w 200000000 -d 1 | grep "RESULT" | tee -a res_pam.txt

# Run Range All
printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./range_test_pam -n $n -r 5 -q $q_all -t 0 -w 200000000 -d 1 | grep "RESULT" | tee -a res_pam.txt

echo
printf "${BLUE}sequential run${NC}\n"
export PARLAY_NUM_THREADS=1
$header ./range_test_pam_seq -n $n -r 5 -q $q_all -t 0 -w 200000000 -d 1 | grep "RESULT" | tee -a res_pam.txt


python3 comp.py
echo
printf "${GREEN}RESULTS (CPAM)${NC}\n"
cat data_cpam.txt
printf "${GREEN}RESULTS (PAM)${NC}\n"
cat data_pam.txt