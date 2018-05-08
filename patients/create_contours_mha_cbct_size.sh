#!/bin/bash

# mask skin from Jihun's data
set -e
for pat in P*; do
    cd $pat
    plastimatch convert --input dicom/ct_extSkinCrop \
                        --output-prefix contours/base_plan_cbctdims \
                        --output-type float
    for s in contours/base_plan_cbctdims/*.mha; do
        plastimatch resample --input "${s}" \
                             --output "${s/.mha/_2.mha}" \
                             --fixed cbct_fraction_1.mha
        plastimatch convert --input "${s/.mha/_2.mha}" \
                            --output-img "${s}" \
                            --output-type uchar
        rm "${s/.mha/_2.mha}"
    done
    cd ..
done