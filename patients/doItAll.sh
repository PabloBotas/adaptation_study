#!/bin/bash

# This is an attempt to show how these scripts should be executed
# The parameters in them, specifically the patients to be processed,
# should be adjusted

###############################################################################
# This step applies only to P02 (for the moment)
# Jihun's data and MCAuto's are misaligned, so we have to fix that
./adjustP02alignment.sh

###############################################################################
# First we can start using the VFs to warp the contours:
./warpMCAutoContours.sh

###############################################################################
# Then we could verify that the CBCTs are in float format
# This is necessary then to transform to ctvolume.dat because the code in
# charge of that is buggy (my fault)
# To see format, run:
plastimatch header file.mha
# Convert to floats with
./convertCBCTtofloat.sh

###############################################################################
# Set HU outside the patient/skin to -1000. -1024 is also fine, gPMC
# sets HU<-1000 to -1000
./maskSkin.sh


###############################################################################
# Finally we can create (or re-create) the dicom directories
./re-createDicoms.sh

