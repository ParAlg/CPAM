#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

export THREADS=`nproc --all`

graph=../inputs/soc-LiveJournal1.adj

header="numactl -i all"
numa_exists="$(numactl --show 2>&1)"
if [[ $numa_exists == *"No NUMA"* ]]; then
    header=""
fi

printf "${GREEN}test simultaneous_updates_queries (CPAM)${NC}\n"
printf "(Results in Figure 11)\n"
make -j $THREADS

printf "${BLUE}Running concurrent queries and updates${NC}\n"
$header ./run_simultaneous_updates_queries-CPAM-CPAM-Diff -query-file query_times_together.csv -update-file update_times_together.csv -s -m -f $graph | grep "RESULT" > results_together.txt
cat results_together.txt

printf "${BLUE}Running only queries${NC}\n"
$header ./run_simultaneous_updates_queries-CPAM-CPAM-Diff -query_only -query-file query_solo.csv -s -m -f $graph | grep "RESULT" > results_query.txt
cat results_query.txt

printf "${BLUE}Running only updates${NC}\n"
$header ./run_simultaneous_updates_queries-CPAM-CPAM-Diff -update_only -update-file update_solo.csv -s -m -f $graph | grep "RESULT" > results_update.txt
cat results_update.txt

printf "${BLUE}Times logged to query_times_together.csv, update_times_together.csv, query_solo.csv, and update_solo.csv, which can be plotted as time series to generate Fig 11. The average throughput and query time are shown above, and in the result files.${NC}\n"
