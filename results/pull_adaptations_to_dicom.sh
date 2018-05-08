#!/bin/bash

set -e

function simply {
    full=${1}
    del=${2}
    full=${full/${del}\//}
    echo ${full}
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
    for d in ${f}/cbct_?; do
        case=$(basename ${d})
        echo "    Case: ${case}"

        cp -r ${ct_dir}/${patient}/dicom/${case/_/_fraction_}_ctdims ${outdir}/${patient}/${case}/cbct
        for adapt_dir in ${d}/adapt_geometric_free ${d}/adapt_geometric_iso_shift; do
            adapt=$(basename ${adapt_dir})
            echo "    Adapt case: ${adapt}"

            for dose in ${adapt_dir}/ast_*.mhd; do
                out="${outdir}/${patient}/${case}"
                echo "        Converting dose: $(simply ${dose} ${adaptive_project})"
                plastimatch convert \
                    --input-dose-img ${dose} \
                    --output-dicom ${out} \
                    --output-type "float" \
                    --referenced-ct ${outdir}/${patient}/${case}/cbct \
                    --metadata "0008,1030"="${adapt}" \
                    --metadata "0008,103e"="${adapt}" \
                    --metadata "300a,00c0"="${adapt}" \
                    --metadata "300a,00c2"="${adapt}" \
                    --metadata "300a,00c3"="${adapt}" > /dev/null
                dose_name=$(basename ${dose})
                dose_name=${dose_name/ast_/}
                dose_name=${dose_name/.mhd/}
                pat_num=${patient%_*}
                mv $(find ${out} -name "*.dcm" -type f | xargs ls -tr | tail -1) "${out}/dose_${pat_num}_${case}_${adapt}_${dose_name}.dcm"
            done
        done
    done
    echo
done
