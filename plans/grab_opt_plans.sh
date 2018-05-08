#!/bin/bash

### This script is intended to grab files from a predefined directory
### The directory is an sshfs volume

origin="/home/gmoc/Desktop/erisroot/users/pb684/adaptive/opt4d_plans"

cd Opt4D
for pat in P04* P08* P10* P14* P16*; do
    mkdir -p ${pat}/optfiles 
    cp $(ls -trd ${origin}/${pat}/job* | tail -n2) ${pat}/optfiles 
    cp ${origin}/${pat}/planfile.pln ${pat}/optfiles
    cp ${origin}/${pat}/results/beam_* ${pat}/optfiles
    cp ${origin}/${pat}/results/DVH.dat ${pat}/optfiles
    cp ${origin}/${pat}/results/dose.dat ${pat}/optfiles
done
cd -
