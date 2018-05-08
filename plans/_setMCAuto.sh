#!/bin/bash

THISPWD=${PWD}
ASTROIDDATA="${THISPWD}/astroid_data/"
OUTDIR="${THISPWD}/Opt4D"

if [ $# -eq 0 ]; then
    input_list=$(ls -d ${OUTDIR}/P10* ${OUTDIR}/P14* ${OUTDIR}/P16*)
else
    input_list="$@"
fi

for out in ${input_list}; do
    echo 
    echo "Setting directory ${out}"
    echo
    pat="$(basename ${out})"
    temp=${pat#*_}
    patid="${temp%${temp#????????????}}"
    temp=${temp/${patid}/}
    file=$(ls -tr ${ASTROIDDATA}/${patid}_*.set | grep "${temp}" | tail -n1)
    astroidID=${file##*__}
    astroidID="${astroidID//[^0-9]/}"
    echo "${patid} ${temp} ${astroidID}"

    # MCAuto files setter and launcher
    cp utils/run_mcauto_template.sh ${out}/run_mcauto.sh
    chmod u+x ${out}/run_mcauto.sh
    sed -i "s,SETID,${astroidID},g" ${out}/run_mcauto.sh
    cd ${out}
    ./run_mcauto.sh
    cd - > /dev/null

    mkdir -p ${out}/astroid_plan
    cp ${ASTROIDDATA}/*${astroidID}*.dose ${out}/astroid_plan
done
