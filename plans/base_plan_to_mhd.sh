#!/bin/bash

## This file runs the passed input patients

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P08*)
else
    input_list="$@"
fi

THISPWD=${PWD}

for dir in ${input_list}; do
    pat=$(basename ${dir})
    pat_name=${pat#*_}
    echo "Running patient ${dir} ..."
    cd ${dir}
    str=$(plastimatch header ${THISPWD}/../patients/${pat}/ct.mha)
    offset=${str#*Origin}
    offset=${offset#*= }
    offset=${offset%Size*}
    n=${str#*Size}
    n=${n#*= }
    n=${n%Spacing*}
    nx=${n%% *}
    nz=${n##* }
    ny=${n/${nx} /}
    ny=${ny/ ${nz}/}
    d=${str#*Spacing}
    d=${d#*= }
    d=${d%Direction*}
    for f in base*/out_*/ast_total.dose; do
        mhd=${f%.*}.mhd
        echo "1, 1, 1, 1, 1, 1, ${nx}, ${ny}, ${nz}" > ${f}.geometry
        echo ${f}
        plastimatch convert \
            --input-dose-ast ${f} \
            --output-dose-img ${mhd}
        sed -i "/Offset.*/c\Offset = ${offset}" ${mhd}
        sed -i "/ElementSpacing.*/c\ElementSpacing = ${d}" ${mhd}
    done
    cd ${THISPWD}
done
