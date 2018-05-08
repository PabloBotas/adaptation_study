#!/bin/bash

set -e

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P* | grep -E "P04|P08|P10|P14|P16")
else
    input_list="$@"
fi

UTILSDIR=${PWD}/utils
for pat in ${input_list}; do
    echo "Adjusting ${pat}"
    patdir=$(basename ${pat})
    patdir=${patdir#*_}
    modifier=${patdir#*_}
    patname=${patdir/$modifier/}
    echo "${patid}"
    tramps=$(ls astroid_data/${patname}*.tramp | grep ${modifier})
    rm -f ${pat}/base_plan/input/tramps/*tramp*
    cp ${tramps} ${pat}/base_plan/input/tramps
    cd ${pat}
    i=1
    for f in base_plan/input/tramps/*.tramp; do
        echo "    - $(basename ${f}) <- beam_${i}.bwf"
        python ${UTILSDIR}/adjustTrampWeights.py \
                --tramp ${f} \
                --weights optfiles/beam_${i}.bwf \
                --fractions 30
        mv ${f} ${f}_modified
        rm ${f}_backup
        ((i++))
    done

    shopt -s nullglob
    for f in cbct_?; do
        rm -f ${f}/input/tramps/*tramp*
        cp base_plan/input/tramps/*.tramp_modified ${f}/input/tramps
    done
    cd - > /dev/null
done

