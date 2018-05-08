#!/bin/bash

# mask skin from Jihun's data
set -e
for pat in P15*; do
    cd $pat
    for c in contours/*fraction*.dcm; do
        num="${c//[!0-9]/}"
        plastimatch convert \
                    --input ${c} \
                    --output-prefix temp_structs/fraction_${num}

        plastimatch convert \
                    --input temp_structs/fraction_${num}/patient.mha \
                    --output-img temp_structs/fraction_${num}/patient.mha \
                    --fixed cbct_fraction_${num}.mha

        mkdir -p temp_mha
        plastimatch mask \
                    --input cbct_fraction_${num}.mha \
                    --mask temp_structs/fraction_${num}/patient.mha \
                    --mask-value -1000 \
                    --output temp_mha/cbct_fraction_${num}.mha

        # if [[ ${pat} != P02* ]]; then
            ### replace cbcts
            mv temp_mha/*.mha .
            rm -rf temp_structs temp_mha
        # fi

    done
    cd ..
done