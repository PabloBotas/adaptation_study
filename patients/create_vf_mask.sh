#!/bin/bash

set -e

function create_files {
    ### This function is only created for convenience to parallelize the
    ### calculation in the background keeping data consistency

    t=${1}
    dmap=${t/.mha/_dmap.mha}
    plastimatch dmap --input ${t} \
                     --output ${dmap}
    # Threshold distance map to include values <= 50 mm and above 0
    rim=${t/.mha/_0.5cm_rim.mha}
    plastimatch threshold --input ${dmap} \
                          --output ${rim} \
                          --range "0.001,5"
    plastimatch multiply ${rim} \
                         ${f}/patient.mha \
                         --output ${rim/.mha/_patient_masked.mha}
    rim=${t/.mha/_1.0cm_rim.mha}
    plastimatch threshold --input ${dmap} \
                          --output ${rim} \
                          --range "0.001,10"
    plastimatch multiply ${rim} \
                         ${f}/patient.mha \
                         --output ${rim/.mha/_patient_masked.mha}
    rim=${t/.mha/_1.5cm_rim.mha}
    plastimatch threshold --input ${dmap} \
                          --output ${rim} \
                          --range "0.001,15"
    rim=${t/.mha/_2.0cm_rim.mha}
    plastimatch threshold --input ${dmap} \
                          --output ${rim} \
                          --range "0.001,20"
    rim=${t/.mha/_3.0cm_rim.mha}
    plastimatch threshold --input ${dmap} \
                          --output ${rim} \
                          --range "0.001,30"
    rm ${dmap}
}

for pat in $(ls -d P01*); do
    echo $pat
    cd $pat
    for f in contours/cbct_1; do
        for t in $(ls ${f}/target*.mha | grep -v "dmap" | grep -v "cm_rim"); do
            create_files ${t} &
        done
    done
    wait
    cd -
done

