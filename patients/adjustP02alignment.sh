#!/bin/bash

set -e

cd P02*
mkdir -p temp_mha
for f in *.mha; do
    plastimatch convert \
                --input "${f}" \
                --xf "transforms/rigid_plas.txt" \
                --output-img "temp_mha/${f}" \
                --fixed "${f}"
done

cd -