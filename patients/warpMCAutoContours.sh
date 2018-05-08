#!/bin/bash

set -e

TEMP_DIR=~/temp/warp_contours.$$.$RANDOM

function empty_temp {
    rm -rf ${TEMP_DIR}/*
}

function clean_up {
    rm -rf ${TEMP_DIR}
    exit $1
}

trap clean_up SIGHUP SIGINT SIGTERM

mkdir -p ${TEMP_DIR}

PLANSDIR=$(realpath "../plans/Opt4D")
for pat in P15*; do
    cd ${pat}
    cp "${PLANSDIR}/${pat}/base_plan/input/ct/rts.dcm" "contours/strs_plan.dcm"
    # for t in transforms/not*xform*.txt; do ## PAT 2
    for t in transforms/week?/xform*.txt; do
        num="${t//[!0-9]/}"
        plastimatch convert \
            --input "contours/strs_plan.dcm" \
            --output-dicom "${TEMP_DIR}" \
            --metadata "0008,103e"="fraction_${num}" \
            --metadata "300a,00c0"="fraction_${num}" \
            --metadata "300a,00c2"="fraction_${num}" \
            --metadata "300a,00c3"="fraction_${num}" \
            --xf "${t}"
        mv ${TEMP_DIR}/*.dcm contours/strs_fraction_${num}.dcm
        empty_temp
    done
    cd -
done

clean_up 0
