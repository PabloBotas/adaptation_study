#!/bin/bash

patient=${1}
cbct=${2}
adapt=${3}

files=$(ls -d Opt4D/${patient}/${cbct}/${adapt}/mask_map.dat \
              Opt4D/${patient}/${cbct}/${adapt}/DoseFrac* \
              Opt4D/${patient}/${cbct}/${adapt}/total.dose \
              Opt4D/${patient}/${cbct}/${adapt}/reopt/*traces.* \
              Opt4D/${patient}/${cbct}/${adapt}/CT.bin)

# Assess
du -ch ${files}

# Clean up
rm ${files}
