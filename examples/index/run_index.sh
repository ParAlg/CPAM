#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

export THREADS=`nproc --all`

num_queries=100000
filename="wiki_small.txt"

header="numactl -i all"
numa_exists="$(numactl --show 2>&1)"
if [[ $numa_exists == *"No NUMA"* ]]; then
    header=""
fi

printf "${GREEN}test inverted index (CPAM)${NC}\n"
printf "(Results in Table 3, row 1)\n"
make -j $THREADS

printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./index -r 5 -q $num_queries -f $filename | grep "RESULT" | tee res.txt

echo
printf "${BLUE}sequential run${NC}\n"
./index_seq -r 5 -q $num_queries -f $filename | grep "RESULT" | tee -a res.txt


printf "${GREEN}test diff-encoded inverted index (CPAM)${NC}\n"
printf "(Results in Table 3, row 2)\n"

printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./index_de -r 5 -q $num_queries -f $filename | grep "RESULT" | tee res_de.txt

echo
printf "${BLUE}sequential run${NC}\n"
./index_de_seq -r 5 -q $num_queries -f $filename | grep "RESULT" | tee -a res_de.txt


echo
printf "${GREEN}test range tree (PAM)${NC}\n"
printf "(Results in Table 3, row 3)\n"

printf "${BLUE}parallel run${NC}\n"
export PARLAY_NUM_THREADS=$THREADS
$header ./index_pam -r 5 -q $num_queries -f $filename | grep "RESULT" | tee res_pam.txt

echo
printf "${BLUE}sequential run${NC}\n"
./index_pam_seq -r 5 -q $num_queries -f $filename | grep "RESULT" | tee -a res_pam.txt


python3 comp.py
echo
printf "${GREEN}RESULTS (CPAM)${NC}\n"
cat data_cpam.txt
printf "${GREEN}RESULTS (CPAM DE)${NC}\n"
cat data_cpam_de.txt
printf "${GREEN}RESULTS (PAM)${NC}\n"
cat data_pam.txt