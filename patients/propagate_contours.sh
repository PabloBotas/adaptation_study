#!/bin/bash

set -e

function create_files {
    input_file=${1}
    output_file=${2}
    vf=${3}
    temp_root=$(dirname ${output_file})
    plastimatch convert --input ${input_file} \
                        --output-img ${output_file} \
                        --output-type "float"
    plastimatch convert --input ${output_file} \
                        --output-img ${output_file} \
                        --xf ${vf} \
                        --output-type "float"
    plastimatch threshold --input ${output_file} \
                          --output ${output_file} \
                          --above 0.5
}

for pat in $(ls -d P01*); do
    cd $pat
    for f in contours/cbct_1; do
        echo ${pat} ${f}
        i=$(basename ${f})
        i=${i/cbct_/}
        # for t in target_Neck_L.mha target_Neck_R.mha oars.mha patient.mha; do
        for t in target.mha oars.mha patient.mha; do
        # for t in oars.mha; do
            input_file=contours/base_plan/${t}
            output_file=${f}/${t}
            vf=transforms/xform_deform_cCBCT${i}-pCT_ctdims.mha
            
            create_files ${input_file} ${output_file} ${vf} &
        done
    done
    wait
    cd -
done

