#!/bin/bash

set -e

input=gpu.in
inputname="${input%.*}"

spotfactor=SPOTFACTOR
fractions=FRACTIONS

/opt/gpromptmc/build/gpromptmc < ${input} > ${inputname}.log

outdir=$(grep "OutputFilePrefix" ${input} | awk -F "/" '{print $NF}')
files=$(ls ${outdir}/DoseAve* | grep -v ast)
/opt/gpmc-tools/1.0/sumDoses ${outdir}/total.dose ${files} >> ${inputname}.log

nx=$(grep "nVoxelsX" ${outdir}/geometry.dat | awk '{print $2}')
ny=$(grep "nVoxelsY" ${outdir}/geometry.dat | awk '{print $2}')
nz=$(grep "nVoxelsZ" ${outdir}/geometry.dat | awk '{print $2}')
for f in $(ls ${outdir}/*.dose | grep -v ast); do
    name=$(basename $f)
    /opt/gpmc-tools/1.0/gpmc2xio ${outdir}/ast_${name} $f $nx $ny $nz $spotfactor $fractions >> ${inputname}.log
done

exit 0