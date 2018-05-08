#!/bin/bash
# set -e

cd Opt4D
INITPWD=${PWD}

for pat in P01* P02* P03* P04*R P05* P07* P08* P10* P14* P15* P16*; do
# for pat in P01*; do
    cd ${pat}
    patnum=${pat%%_*}
    temp=${pat/${patnum}_/}
    patname=${temp%%_*}
    patid=${patnum}_${patname}
    target_case=${pat/${patid}/}
    echo $pat

    for fracdir in cbct_?; do
        cd ${fracdir}
        i=${fracdir##*_}
        frac=cbct_${i}
        echo "  $fracdir"

        for method in "geometric"; do
            for strat in "free" "iso_shift" "range_shifter" "range_shifter iso_shift"; do
            # for strat in "free"; do
                x="${strat// /_}"
                adaptcase="adapt_${method}_${x}"
                log=${adaptcase}/${adaptcase}.log

                mkdir -p "${adaptcase}"
                echo "    $strat"
                contours=../../../../patients/${patid}/contours
                transforms=../../../../patients/${patid}/transforms

                /opt/utils/adaptive-build/adaptive \
                        --plan_dir    ../base_plan \
                        --adapt_dir   . \
                        --cbct        ${fracdir}.mha \
                        --vf_mask     1 ${contours}/base_plan/target${target_case}.mha \
                                      1 ${contours}/base_plan/target${target_case}_1.0cm_rim_patient_masked.mha \
                                      ${contours}/base_plan/target${target_case}_1.5cm_rim.mha \
                                      ${contours}/base_plan/target${target_case}_2.0cm_rim.mha \
                                      ${contours}/base_plan/target${target_case}_3.0cm_rim.mha \
                        --target_mask ${contours}/${frac}/target${target_case}.mha  \
                        --target_rim  ${contours}/${frac}/target${target_case}_1.0cm_rim_patient_masked.mha \
                        --oars        ${contours}/${frac}/oars.mha \
                        --vf          ${transforms}/xform_deform_cCBCT${i}-pC*.mha \
                        --optdir      ${adaptcase}/reopt \
                        --outplan     ${adaptcase} \
                        --out_vf      --out_shifts \
                        --constraint  "${strat}" \
                        --method      "${method}" \
                        --traces_ct   --traces_cbct \
                        --sf_dose 1000000 \
                        --sim_adapt \
                        --plan_dose  ../base_plan/out_base_plan/DosePlan.total.dose \
                        --vf_report  ${adaptcase}_${patnum}_${fracdir}.pdf \
                        > ${log}
                        # > >(tee ${log}) 2> >(tee -a ${log} >&2)

                if [[ $? == 0 ]]; then
                    for f in $(ls -d vf_mask_map.dat \
                               ${adaptcase}/DoseFrac* \
                               ${adaptcase}/ast_DoseFrac*.mhd \
                               ${adaptcase}/ast_DoseFrac*.raw \
                               ${adaptcase}/total.dose \
                               ${adaptcase}/reopt/*traces.* \
                               ${adaptcase}/CT.bin); do
                        if [ -f ${f} ]; then
                            rm -f ${f}
                        fi
                    done
                fi
            done
        done
        cd ..
    done
    cd ${INITPWD}
done
cd ..
