#!/bin/bash


for pat in P0[4,5,7,8]*; do
    cd ${pat}
    # for cbct in cbct_fraction_?.mha; do
    #     plastimatch resample --fixed ${fixed} \
    #                          --input ${cbct} \
    #                          --output ${cbct/.mha/_ctdims.mha} \
    #                          --default-value -1000
    # done
    for vf in transforms/xform*.mha; do
        plastimatch resample --fixed ct.mha --input ${vf} --output ${vf/.mha/_ctdims.mha}
    done
    # cp -r contours contours_original
    # for v in contours/*/*.mha; do
        # plastimatch resample --fixed ${fixed} --input "${v}" --output "${v}"
    # done
    cd -
done

exit 0