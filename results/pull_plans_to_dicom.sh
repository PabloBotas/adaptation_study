#!/bin/bash

set -e

function simply {
    full=${1}
    del=${2}
    full=${full/${del}\//}
    echo ${full}
}

function find_ct_slices {
    dir=${1}
    exclude=${2}
    if [ ! -z "${exclude}" ]; then
        slices=""
        for f in $(ls ${dir}/*); do
            if [[ $(basename ${f}) != *"${exclude}"* ]]; then
                slices="${slices} ${f}"
            fi
        done
    else
        slices=$(ls ${dir}/*)
    fi
    echo ${slices}
}

plan_dir=$(realpath "../plans/Opt4D")
astroid_dir=$(realpath "../plans")
ct_dir=$(realpath "../patients")
outdir=$(realpath "data/dicom")
adaptive_project=$(realpath ..)

echo "Paths printed relative to ${adaptive_project}"
echo

if [ $# -eq 0 ]; then
    input_list=$(ls -d ${plan_dir}/P*)
else
    input_list="$@"
fi

for f in ${input_list}; do
    patient=$(basename ${f})
    patient_name=${patient#*_}
    echo "Working on patient ${patient}"

    # Loop around cases converting and copying data to DICOM
    for d in ${f}/base_plan ${f}/cbct_*; do
        case=$(basename ${d})
        echo "    Case: ${case}"

        if [[ ${case} == "astroid_plan" ]] || [[ ${case} == "base_plan" ]]; then
            outcase="plan"
        else
            outcase=${case}
        fi

        ## Copy CT data
        mkdir -p ${outdir}/${patient}/${outcase}/ct
        if [[ ${case} == "base_plan" ]]; then
            echo "        Copying CT data: $(simply ${plan_dir} ${adaptive_project})/${patient}${case}/input/ct"
            slices=$(find_ct_slices ${plan_dir}/${patient}/${case}/input/ct)
            cp ${slices} ${outdir}/${patient}/${outcase}/ct
        elif [[ ! ${case} == "astroid_plan" ]]; then
            temp=${case/_/_fraction_}
            echo "        Copying CT data: $(simply ${ct_dir} ${adaptive_project})/${patient}/dicom/${temp}"
            slices=$(find_ct_slices ${ct_dir}/${patient}/dicom/${temp} str)
            cp ${slices} ${outdir}/${patient}/${outcase}/ct
            temp=${case/cbct/strs_fraction}
            temp=${temp}.dcm
            echo "        Contours set to: $(simply ${ct_dir} ${adaptive_project})/${patient}/contours/${temp}"
            cp ${ct_dir}/${patient}/contours/${temp} ${outdir}/${patient}/${outcase}/ct
        fi

        ## Copy doses
        if [[ ${case} == "astroid_plan" ]]; then
            dose=$(find ${d} -type f -name '*.dose' | sort | tail -n 1)
        else
            dose="${d}/out_${case}/ast_total.dose"
        fi

        ## Create .geometry file needed by plastimatch
        str=$(head -n1 ${astroid_dir}/astroid_data/${patient_name}*.set)
        str=${str/calc-vol\ \"/}
        str=${str/\"/}
        echo ${str} > ${dose}.geometry

        out="${outdir}/${patient}/${outcase}"
        echo "        Converting dose: $(simply ${dose} ${adaptive_project})"
        plastimatch convert \
            --input-dose-ast ${dose} \
            --output-dicom ${out} \
            --output-type "float" \
            --referenced-ct ${outdir}/ct \
            --metadata "0008,1030"="${case}" \
            --metadata "0008,103e"="${case}" \
            --metadata "300a,00c0"="${case}" \
            --metadata "300a,00c2"="${case}" \
            --metadata "300a,00c3"="${case}" > /dev/null
        mv $(find ${out} -name "*.dcm" -type f | xargs ls -tr | tail -1) "${out}/dose_${patient}_${case}.dcm"
    done
    echo
done
