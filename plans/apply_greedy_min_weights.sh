#!/bin/bash

## This file runs the passed input patients

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P01*)
else
    input_list="$@"
fi

THISPWD=${PWD}

for dir in ${input_list}; do
    pat=$(basename ${dir})
    pat_name=${pat#*_}
    echo "Running patient ${dir} ..."
    cd ${dir}
    for tramp in base_plan/input/tramps/*.tramp; do
        echo "    ${tramp}"
        /opt/astroid_database/tramp_greedy_MediumSpots.sh ${tramp} ${tramp}_modified
    done
    # for t in cbct_?/input/tramps; do
    #     cp base_plan/input/tramps/*.tramp_modified ${t}
    # done
    cd ${THISPWD}
done
