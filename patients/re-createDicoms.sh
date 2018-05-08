#!/bin/bash

set -e

for pat in P*; do
    cd $pat
    for f in *.mha; do
        outname=${f/.mha/}
        rm -rf dicom/${outname}
        plastimatch convert \
                    --input "${f}" \
                    --output-dicom "dicom/${outname}"   \
                    --metadata "0008,1030"="${outname}" \
                    --metadata "0008,103e"="${outname}" \
                    --metadata "300a,00c0"="${outname}" \
                    --metadata "300a,00c2"="${outname}" \
                    --metadata "300a,00c3"="${outname}"

        if [[ ${f} == ct* ]]; then
            cp contours/strs_plan.dcm dicom/${outname}
        else
            cnt=${f/cbct/strs}
            cnt=${cnt/.mha/.dcm}
            cp contours/${cnt} dicom/${outname}
        fi
    done
    cd -
done