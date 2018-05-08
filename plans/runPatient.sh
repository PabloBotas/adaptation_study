#!/bin/bash

## This file runs the passed input patients

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P04*R Opt4D/P10* Opt4D/P08*)
else
    input_list="$@"
fi

THISPWD=${PWD}

shopt -s nullglob
for dir in ${input_list}; do
    pat=$(basename ${dir})
    echo "Running patient ${dir} ..."
    cd ${dir}
    for case in base_plan cbct_?; do
        echo "    ${case}"
        cd ${case}
        ./gpu_run.sh
        cd - > /dev/null
    done
    cd ${THISPWD}
done
