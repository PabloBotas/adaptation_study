#!/bin/bash

## This file runs the passed input patients

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P01* Opt4D/P02* Opt4D/P03* Opt4D/P04*R Opt4D/P05* Opt4D/P07* Opt4D/P10* Opt4D/P14* Opt4D/P15* Opt4D/P16*)
    #input_list=$(ls -d Opt4D/P01*)
else
    input_list="$@"
fi

THISPWD=${PWD}

for dir in Opt4D/P01* Opt4D/P02* Opt4D/P03* Opt4D/P04*R Opt4D/P05* Opt4D/P07* Opt4D/P10* Opt4D/P14* Opt4D/P15* Opt4D/P16*; do
# for dir in ${input_list}; do
    pat=$(basename ${dir})
    pat_name=${pat#*_}
    echo "Running patient ${dir} ..."
    cd ${dir}
    # for cbct in cbct_1 cbct_2 cbct_3 cbct_4 cbct_5 cbct_6; do
    for cbct in cbct_?; do
        echo "    ${cbct}"
        cd ${cbct}
        for method in "cold_spots"; do
        # for method in "geometric" "cold_spots"; do
            for strat in "free" "iso_shift" "range_shifter" "range_shifter iso_shift"; do
            # for strat in "free"; do
                x="${strat// /_}"
                adaptcase="adapt_${method}_${x}"

                echo "        ${adaptcase}"
                str=$(plastimatch header ${cbct}.mha)
                spacing=${str#*Spacing =}
                spacing=${spacing%Direction*}
                offset=${str#*Origin =}
                offset=${offset%Size*}
                n=${str#*Size}
                n=${n#*= }
                n=${n%Spacing*}
                nx=${n%% *}
                nz=${n##* }
                ny=${n/${nx} /}
                ny=${ny/ ${nz}/}
                for f in ${adaptcase}/ast_total.dose; do
                    mhd=${f%.*}.mhd
                    echo "1, 1, 1, 1, 1, 1, ${nx}, ${ny}, ${nz}" > ${f}.geometry
                    # echo ${f}
                    plastimatch convert \
                        --input-dose-ast ${f} \
                        --output-dose-img ${mhd} > /dev/null
                    ret_code=$?
                    if [ $ret_code != 0 ]; then
                        printf "Error : [%d] when converting ${f}" $ret_code
                        exit $ret_code
                    else
                        rm ${f}
                    fi
                    sed -i "/Offset.*/c\Offset = ${offset}" ${mhd}
                    sed -i "/ElementSpacing.*/c\ElementSpacing = ${spacing}" ${mhd}
                done
            done
        done
        cd ..
    done
    cd ${THISPWD}
done
