#!/bin/bash

set -e

for pat in P07* P08* P10* P14* P15* P16*; do
    echo ""
    echo ${pat}
    echo ""
    cd ${pat}/contours
    i=1
    for t in ../transforms/xform*.mha; do
        echo "Step 1"
        plastimatch convert \
            --input-prefix "base_plan" \
            --output-prefix "temp"
        echo "Step 2"
        for f in temp/*.mha; do
            plastimatch convert \
                --input "${f}" \
                --output-img "${f}" \
                --output-type "float"
            plastimatch convert \
                --input "${f}" \
                --output-img "${f}" \
                --xf "${t}"
            plastimatch threshold \
                --input "${f}" \
                --output "${f}" \
                --above "0.5"
        done
        echo "Step 3"
        plastimatch convert \
            --input-prefix "temp" \
            --output-ss-img "structs_ss_cbct_${i}.mha" \
            --output-ss-list "structs_ss_cbct_${i}.txt"
        rm -rf temp
        ((i++))
    done
    cd -
done


# for pat in P08*; do
#     cd ${pat}
#     i=1
#     for t in transforms/xform*.txt; do
#         for f in contours/base_plan/targe*.mha; do
#             plastimatch convert \
#                 --input "${f}" \
#                 --output-img "contours/cbct_${i}/$(basename ${f})" \
#                 --fixed "${f}" \
#                 --output-type "uchar" \
#                 --xf "${t}"
#         done &
#         ((i++))
#     done
#     cd -
# done

