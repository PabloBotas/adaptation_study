#!/bin/bash


cd Opt4D
for pat in $(ls -d P* | grep -v "_L" | grep -v "P08"); do
    echo ${pat}
    ## ADAPTED STRATEGIES =================================================
    for strat in ${pat}/cummulative/adapt*; do
        echo ${strat}
        ## Propagate dose to same set of contours with inverted VF
        i=1
        for f in ${pat}/cbct_?/$(basename ${strat})/ast_total.mhd; do
            plastimatch convert --input ${f} \
                                --output-img ${strat}/warped_cbct_${i}.mhd \
                                --xf ../../patients/${pat/_Neck_R/}/transforms/inverted*${i}*.mha > /dev/null &
            ((i++))
        done
        wait
        plastimatch average --output ${strat}/cummulative_dose.mhd \
                            ${strat}/warped_cbct_?.mhd
        if [[ $? == 0 ]]; then
            rm -f ${strat}/warped_cbct_?.mhd
        fi
    done

    ## UNADAPTED STRATEGIES =================================================
    ## Propagate dose to same set of contours with inverted VF
    # i=1
    # for f in ${pat}/cbct_?/out_cbct_?/ast_total.mhd; do
    #     plastimatch convert --input ${f} \
    #                         --output-img ${pat}/cummulative/warped_cbct_${i}.mhd \
    #                         --xf ../../patients/${pat/_Neck_R/}/transforms/inverted*${i}*.mha &
    #     ((i++))
    # done
    # wait
    # plastimatch average --output ${pat}/cummulative/cummulative_unadapted.mhd \
    #                     ${pat}/cummulative/warped_cbct_?.mhd
    # if [[ $? == 0 ]]; then
    #     rm -f ${pat}/cummulative/warped_cbct_?.mhd
    # fi
    
done
cd ..


