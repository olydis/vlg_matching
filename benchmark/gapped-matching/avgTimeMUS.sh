#!/bin/bash
cat results/all.txt | grep 'TIMING' | awk '{total+=$3} END {print total/NR}'
