#!/bin/bash

## This file runs the passed input patients

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P*)
else
    input_list="$@"
fi

for dir in ${input_list}; do
    cd ${dir}
    mkdir -p optfiles
    mv job.* planfile.pln results/* optfiles 2> /dev/null || true
    rm -rf results
    cd - > /dev/null
done
