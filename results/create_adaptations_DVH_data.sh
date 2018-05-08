#!/bin/bash

set -e

cd data

LOGDIR=$(realpath logs)
DVHDIR=$(realpath dvh)
PLANSDIR=$(realpath ../../plans/Opt4D)
STRDIR=$(realpath ../../patients)
TEMPDIR=$(realpath temp)
mkdir -p ${LOGDIR}
mkdir -p ${DVHDIR}
mkdir -p ${TEMPDIR}

if [ $# -eq 0 ]; then
    input_list=$(ls -d ${PLANSDIR}/P* | grep -v "_L" | grep -v "P08")
    # input_list=$(ls -d ${PLANSDIR}/P* | grep -E "P10|P16")
else
    input_list="$@"
fi

function simply {
    full=${1}
    del=${2}
    full=${full/${del}\//}
    echo ${full}
}

for patpath in ${input_list}; do
    pat=$(basename ${patpath})
    echo "Creating DVH for ${pat}"
    cd ${patpath}
    for cbct in cbct_?; do
        for method in "cold_spots" "geometric"; do
        # for method in "cold_spots"; do
            for strat in "free" "iso_shift" "range_shifter" "range_shifter iso_shift"; do
            # for strat in "free"; do
                strat="${strat// /_}"
                adaptcase="adapt_${method}_${strat}"

                ### ADAPTATIONS
                fraccase="${cbct}_${adaptcase}"
                rm -f ${LOGDIR}/${pat}_${fraccase}.log
                for dose in ${cbct}/${adaptcase}/ast_total.mhd; do
                    echo "    -${fraccase}: $(simply ${dose} ${PLANSDIR})"
                    # plastimatch stats "${dose}"
                    plastimatch dvh \
                        --input-dose    "${dose}" \
                        --input-ss-img  "${STRDIR}/${pat/_Neck_R/}/contours/structs_ss_${cbct}.mha" \
                        --input-ss-list "${STRDIR}/${pat/_Neck_R/}/contours/structs_ss_${cbct}.txt" \
                        --output-csv    "${DVHDIR}/${pat}_${fraccase}.dvh" \
                        --dose-units    "gy"                               \
                        --num-bins      "901"                              \
                        --bin-width     "0.1"                              \
                        --cumulative >> ${LOGDIR}/${pat}_${fraccase}.log &
                done

                ### CUMULATIVE ADAPTATIONS
                cumulative="${patpath}/cumulative"
                dose="${cumulative}/${adaptcase}/cummulative_dose.mhd"
                fraccase="cumulative_${adaptcase}"
                plastimatch dvh \
                       --input-dose    "${dose}" \
                       --input-ss-img  "${STRDIR}/${pat/_Neck_R/}/contours/structs_ss_base_plan.mha" \
                       --input-ss-list "${STRDIR}/${pat/_Neck_R/}/contours/structs_ss_base_plan.txt" \
                       --output-csv    "${DVHDIR}/${pat}_${fraccase}.dvh" \
                       --dose-units    "gy"                               \
                       --num-bins      "901"                              \
                       --bin-width     "0.1"                              \
                       --cumulative >> ${LOGDIR}/${pat}_${fraccase}.log &

                ### CUMULATIVE NON-ADAPTATIONS
                cumulative="${patpath}/cumulative"
                dose="${cumulative}/cummulative_unadapted.mhd"
                fraccase="cumulative"
                plastimatch dvh \
                       --input-dose    "${dose}" \
                       --input-ss-img  "${STRDIR}/${pat/_Neck_R/}/contours/structs_ss_base_plan.mha" \
                       --input-ss-list "${STRDIR}/${pat/_Neck_R/}/contours/structs_ss_base_plan.txt" \
                       --output-csv    "${DVHDIR}/${pat}_${fraccase}.dvh" \
                       --dose-units    "gy"                               \
                       --num-bins      "901"                              \
                       --bin-width     "0.1"                              \
                       --cumulative >> ${LOGDIR}/${pat}_${fraccase}.log &
            done
            wait
        done
    done
    cd -
done

cd ..
