#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

cd ../other_systems/aspen/
make run_static_algorithm
cd ../../run_static_algorithms

printf "${BLUE}Running CPAM algorithms.${NC}\n"
python3 run_cpam.py
printf "${BLUE}Running Aspen algorithms.${NC}\n"
python3 run_aspen.py
printf "${BLUE}Running eval.${NC}\n"
python3 eval.py
