#!/bin/bash
set -e

cd Opt4D
# for method in "gpmc_dij"; do
for method in "geometric"; do
    # for strat in "free" "virt_range_shifter" "iso_shift" "virt_range_shifter iso_shift"; do
    # for strat in "free" "iso_shift" "range_shifter" "range_shifter iso_shift"; do
    # for strat in "free" "iso_shift" "range_shifter"; do
    for strat in $1; do
        x="${strat// /_}"
        adaptcase="adapt_${method}_${x}"
        log=${adaptcase}/${adaptcase}.log
        INITPWD=${PWD}
        for pat in P15*; do
            cd ${pat}
            patnum=${pat%%_*}
            temp=${pat/${patnum}_/}
            patname=${temp%%_*}
            patid=${patnum}_${patname}
            target_case=${pat/${patid}/}
            echo $pat
            for fracdir in cbct_?; do
            # for fracdir in warped_ct_1 warped_ct_3 warped_ct_6; do
            # for fracdir in cbct_6; do
                cd ${fracdir}
                i=${fracdir##*_}
                frac=cbct_${i}
                mkdir -p "${adaptcase}"
                echo "    $frac"
                echo "        $strat"
                contours=../../../../patients/${patid}/contours
                transforms=../../../../patients/${patid}/transforms
                /opt/utils/adaptive-build/adaptive \
                        --plan_dir   ../base_plan \
                        --adapt_dir  . \
                        --cbct       ${fracdir}.mha \
                        --vf_mask     1 ${contours}/base_plan/target${target_case}.mha \
                                      1 ${contours}/base_plan/target${target_case}_1.0cm_rim_patient_masked.mha \
                                      ${contours}/base_plan/target${target_case}_1.5cm_rim.mha \
                                      ${contours}/base_plan/target${target_case}_2.0cm_rim.mha \
                                      ${contours}/base_plan/target${target_case}_3.0cm_rim.mha \
                        --target_mask ${contours}/${frac}/target${target_case}.mha  \
                        --target_rim  ${contours}/${frac}/target${target_case}_1.0cm_rim_patient_masked.mha \
                        --oars        ${contours}/${frac}/oars.mha \
                        --vf         ${transforms}/xform_deform_cCBCT${i}-pCT_ctdims.mha \
                        --optdir     ${adaptcase}/reopt \
                        --outplan    ${adaptcase} \
                        --out_vf     --out_shifts \
                        --constraint "${strat}" \
                        --method     "${method}" \
                        --traces_ct  --traces_cbct \
                        --sf_dose 1000000 \
                        --sim_adapt  \
                        --plan_dose  ../base_plan/out_base_plan/DosePlan.total.dose \
                        --vf_report  ${adaptcase}_${patnum}_${fracdir}.pdf \
                        > ${log}
                        # > >(tee ${log}) 2> >(tee -a ${log} >&2)
                cd ..
            done
            cd ${INITPWD}
        done    
    done
done
cd ..
