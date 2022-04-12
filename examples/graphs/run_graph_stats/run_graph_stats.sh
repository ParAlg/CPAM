#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color


printf "${BLUE}Building GBBS Inputs.${NC}\n"

gbbs_converter="../other_systems/gbbs/utils/converter"

echo "Building converter (.adj -> .bytepda)"
cd ../other_systems/gbbs/utils
make converter
cd ../../../run_graph_stats

cd ../inputs
for fname in *.adj
do
  compressed_fname=${fname%.adj}.bytepda
  if ! [ -e $compressed_fname ]
  then
    $gbbs_converter -s -bs 128 -rounds 1 -m -enc bytepd-amortized -o $compressed_fname $fname
  fi
done

cd ../run_graph_stats


printf "${BLUE}Building Aspen's memory_footprint tool.${NC}\n"
cd ../other_systems/aspen

make memory_footprint

memory_footprint=../other_systems/aspen/memory_footprint

cd ../../run_graph_stats

printf "${BLUE}Building CPAM-based and PAM-based graph representations.${NC}\n"
make -j $THREADS

printf "${BLUE}Running stats computations.${NC}\n"

for fname in ../inputs/*.adj
do
  printf "${BLUE}Running on ${fname}.${NC}\n"
  name="${fname##*/}"
  output=${name%.adj}.output
  ./GraphStatsPAM-PAM -rounds 1 -s -m $fname | grep "csv" | tee $output > /dev/null
  ./GraphStatsPAM-CPAM -rounds 1 -s -m $fname | grep "csv" | tee -a $output > /dev/null
  ./GraphStatsCPAM-CPAM -rounds 1 -s -m $fname | grep "csv" | tee -a $output > /dev/null
  ./GraphStatsCPAM-CPAM-Diff -rounds 1 -s -m $fname | grep "csv" | tee -a $output > /dev/null
  $memory_footprint -rounds 1 -s -m -f $fname | grep "calculated" | tee -a $output > /dev/null

  du -b ${fname%.adj}.bytepda | tee -a $output > /dev/null
done

python3 run_eval.py
