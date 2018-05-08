#!/bin/bash

set -e

TEMP_DIR=~/temp/$$.$RANDOM
TEMP_FILE=${TEMP_DIR}/cbct_temp.mha

function empty_temp {
    rm -rf ${TEMP_DIR}/*
}

function clean_up {
    rm -rf ${TEMP_DIR}
    exit $1
}

trap clean_up SIGHUP SIGINT SIGTERM

mkdir -p ${TEMP_DIR}

for f in P*; do
    cd "${f}"
    i=1
    for s in cbct_fraction*; do
        plastimatch convert \
            --input "${s}" \
            --output-img "${TEMP_FILE}" \
            --output-type "float"
        mv ${TEMP_FILE} ${s}
        empty_temp
    done
    cd ..
done

clean_up 0
