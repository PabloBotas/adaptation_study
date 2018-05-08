#!/bin/bash

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P* | grep -E "P10|P14|P16")
else
    input_list="$@"
fi

THISPWD=${PWD}
PATDIR=${THISPWD}/../patients

for dir in ${input_list}; do
    pat=$(basename ${dir})
    patnum=${pat%%_*}
    temp=${pat#*_}
    patname=${temp%%_*}
    pat=${patnum}_${patname}
    echo $pat
    i=1
    for f in ${PATDIR}/${pat}/cbct_*ctdims.mha; do
        mkdir -p ${dir}/cbct_${i}/input
        cp -r $(ls -d ${dir}/base_plan/input/* | grep -v cradle | grep -vw ct) ${dir}/cbct_${i}/input
        cp ${f} ${dir}/cbct_${i}/cbct_${i}.mha
        mha2ctvolume ${f} ${dir}/cbct_${i}/input/ctbinary/ctvolume.dat
        ((i++))
    done
done
