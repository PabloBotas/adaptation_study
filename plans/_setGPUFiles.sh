#!/bin/bash

set -e

spotfactor=1000000
fractions=30

if [ $# -eq 0 ]; then
    input_list=$(ls -d Opt4D/P0[4,8]* Opt4D/P1[0,4,6]*)
else
    input_list="$@"
fi

shopt -s nullglob
for dir in ${input_list}; do
    pat=$(basename ${dir})
    for f in ${dir}/base_plan ${dir}/cbct_?; do
        plan=$(basename ${f})
        # GPU launcher
        cp utils/run_gpu_template.sh ${f}/gpu_run.sh
        chmod u+x ${f}/gpu_run.sh
        sed -i "s,FRACTIONS,${fractions},g" ${f}/gpu_run.sh
        sed -i "s,SPOTFACTOR,${spotfactor},g" ${f}/gpu_run.sh
        # GPU plan files setter
        cp utils/gpu_template.in ${f}/gpu.in
        sed -i "s,PATIENTDIR,.,g" ${f}/gpu.in
        sed -i "s,CTVOLUME_FILE,input/ctbinary/ctvolume.dat,g" ${f}/gpu.in
        sed -i "s,OUTPUT,out_${plan},g" ${f}/gpu.in
        sed -i "s,SPOTFACTOR,${spotfactor},g" ${f}/gpu.in
    done
done

# #### Change MCAuto files to adapt to new CT dimensions
# function set_field {
#     field=${1}
#     value=${2}
#     file=${3}
#     sed -i "s,^${field}.*,${field}= ${value}," $file
# }
# function get_field {
#     field=${1}
#     file=${2}
#     line=$(grep "${field}" ${file})
#     line="${line/\ mm/}"
#     value="${line/${field}=\ /}"
#     echo "${value}"
# }
# function get_vectorfield {
#     field=${1}
#     file=${2}
#     line=$(grep "${field}" ${file})
#     line="${line/\ mm/}"
#     value="${line/${field}=\ /}"
#     echo ${value} | awk '{print $2}'
# } 

# cbct_cases=${input_list}
# for dir in ${cbct_cases}; do
#     pat=$(basename ${dir})

#     case=${dir}/cbct_1
#     echo ${case}
#     cd ${case}
#     for beamstr in $(ls -d input/* | grep -v ct | grep -v tramps | grep -v cradle); do
#         beam=$(basename ${beamstr})
#         cbct=$(basename ${case}).mha
#         echo ${beam}
#         original=${beamstr}/run/MCAUTO_DICOM.txt
#         if [ -f ${original/.txt/_backup.txt} ]; then
#             echo "File already adapted, skipping :)"
#             continue
#         fi
#         cp ${original} ${original/.txt/_backup.txt}
#         new=${beamstr}/run/MCAUTO_DICOM.txt

#         oImgCtrX=$(get_field "d:Rt/CT/ImgCenterX"                    ${original})
#         oImgCtrY=$(get_field "d:Rt/CT/ImgCenterY"                    ${original})
#         oImgCtrZ=$(get_field "d:Rt/CT/ImgCenterZ"                    ${original})
#         oSpacingX=$(get_field "d:Rt/CT/PixelSpacing0"                ${original})
#         oSpacingY=$(get_field "d:Rt/CT/PixelSpacing1"                ${original})
#         oSpacingZ=$(get_vectorfield "dv:Rt/CT/SliceThicknessSpacing" ${original})
#         oSizeX=$(get_field "i:Rt/CT/Rows"                            ${original})
#         oSizeY=$(get_field "i:Rt/CT/Columns"                         ${original})
#         oSizeZ=$(get_vectorfield "uv:Rt/CT/SliceThicknessSections"   ${original})
#         oDose2CTX=$(get_field "d:Rt/dose/DoseGridCenter2CTOriginX"   ${original})
#         oDose2CTY=$(get_field "d:Rt/dose/DoseGridCenter2CTOriginY"   ${original})
#         oDose2CTZ=$(get_field "d:Rt/dose/DoseGridCenter2CTOriginZ"   ${original})

#         oOriginX=$(echo "scale=3; ${oImgCtrX} + 1/2*(1-${oSizeX})*${oSpacingX}" | bc | awk '{printf "%f", $0}')
#         oOriginY=$(echo "scale=3; ${oImgCtrY} + 1/2*(1-${oSizeY})*${oSpacingY}" | bc | awk '{printf "%f", $0}')
#         oOriginZ=$(echo "scale=3; ${oImgCtrZ} + 1/2*(1-${oSizeZ})*${oSpacingZ}" | bc | awk '{printf "%f", $0}')

#         header=$(plastimatch header ${cbct})
#         origin=$(echo "${header}" | grep "Origin")
#         originX=$(echo ${origin} | awk '{print $3}')
#         originY=$(echo ${origin} | awk '{print $4}')
#         originZ=$(echo ${origin} | awk '{print $5}')
#         size=$(echo "${header}" | grep "Size")
#         sizeX=$(echo ${size} | awk '{print $3}')
#         sizeY=$(echo ${size} | awk '{print $4}')
#         sizeZ=$(echo ${size} | awk '{print $5}')
#         spacing=$(echo "${header}" | grep "Spacing")
#         spacingX=$(echo ${spacing} | awk '{print $3}')
#         spacingY=$(echo ${spacing} | awk '{print $4}')
#         spacingZ=$(echo ${spacing} | awk '{print $5}')

#         imgcenterX=$(echo "scale=3; ${originX} - ${spacingX}/2 + ${spacingX}*${sizeX}/2" | bc | awk '{printf "%f", $0}')
#         imgcenterY=$(echo "scale=3; ${originY} - ${spacingY}/2 + ${spacingY}*${sizeY}/2" | bc | awk '{printf "%f", $0}')
#         imgcenterZ=$(echo "scale=3; ${originZ} - ${spacingZ}/2 + ${spacingZ}*${sizeZ}/2" | bc | awk '{printf "%f", $0}')
#         dose2CTX=$(echo "scale=3; ${originX} - ${oOriginX} + ${oDose2CTX}" | bc | awk '{printf "%f", $0}')
#         dose2CTY=$(echo "scale=3; ${originY} - ${oOriginY} + ${oDose2CTY}" | bc | awk '{printf "%f", $0}')
#         dose2CTZ=$(echo "scale=3; ${originZ} - ${oOriginZ} + ${oDose2CTZ}" | bc | awk '{printf "%f", $0}')

#         set_field "d:Rt/CT/ImgCenterX" "${imgcenterX} mm" ${new}
#         set_field "d:Rt/CT/ImgCenterY" "${imgcenterY} mm" ${new}
#         set_field "d:Rt/CT/ImgCenterZ" "${imgcenterZ} mm" ${new}
#         set_field "d:Rt/CT/PixelSpacing0"            "${spacingX} mm" ${new}
#         set_field "d:Rt/CT/PixelSpacing1"            "${spacingY} mm" ${new}
#         set_field "dv:Rt/CT/SliceThicknessSpacing" "1 ${spacingZ} mm" ${new}
#         set_field "i:Rt/CT/Columns"                   "${sizeX}" ${new}
#         set_field "i:Rt/CT/Rows"                      "${sizeY}" ${new}
#         set_field "uv:Rt/CT/SliceThicknessSections" "1 ${sizeZ}" ${new}
#     done
#     cd - > /dev/null


#     for case in $(ls -d ${dir}/cbct_? | grep -v cbct_1); do
#         echo ${case}
#         cd ${case}
#         for beamstr in $(ls -d input/* | grep -v ct | grep -v tramps | grep -v cradle); do
#             cp ../cbct_1/${beamstr}/run/MCAUTO_DICOM.txt ${beamstr}/run/MCAUTO_DICOM.txt
#         done
#         cd - > /dev/null 
#     done

# done
