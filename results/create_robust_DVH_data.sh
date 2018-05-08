#!/bin/bash

set -e

cd data

STRDIR=$(realpath ../../patients)
LOGDIR="logs"
DVHDIR="dvh"
PLANSDIR=$(realpath ../../plans/Opt4D)
TEMPDIR="temp"
mkdir -p ${STRDIR}
mkdir -p ${LOGDIR}
mkdir -p ${DVHDIR}
mkdir -p ${TEMPDIR}

function simply {
    full=${1}
    del=${2}
    full=${full/${del}\//}
    echo ${full}
}

for patpath in $(ls -d ${PLANSDIR}/P* | grep -E "P01"); do
    patient=$(basename ${patpath})
    echo "Creating DVH for ${patient}"
    for dose in ${patpath}/*/out_robust_cbct*/ast_total.mhd; do
        case=$(dirname ${dose})
        case=$(dirname ${case})
        case=$(basename ${case})
        echo "    -${case}: $(simply ${dose} ${PLANSDIR})"
        plastimatch dvh \
            --input-dose    "${dose}"     \
            --input-ss-img  "${STRDIR}/${patient}/contours/structs_ss_${case}.mha" \
            --input-ss-list "${STRDIR}/${patient}/contours/structs_ss_${case}.txt" \
            --output-csv    "${DVHDIR}/${patient}_${case}_robust.dvh"           \
            --dose-units    "gy"                                         \
            --num-bins      "901"                                        \
            --bin-width     "0.1"                                        \
            --cumulative >> ${LOGDIR}/${patient}_${case}.log
    done
done

cd ..
