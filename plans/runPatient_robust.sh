#!/bin/bash

## This file runs the passed input patients

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P01*)
else
    input_list="$@"
fi

THISPWD=${PWD}

shopt -s nullglob
for dir in ${input_list}; do
    pat=$(basename ${dir})
    echo "Running patient ${dir} ..."
    cd ${dir}
    for i in 1 2 3 4 5 6; do
        case=cbct_${i}
        echo "    ${case}"
        cd ${case}
        ./robust_gpu_run.sh
        cd - > /dev/null
    done
    cd ${THISPWD}
done
