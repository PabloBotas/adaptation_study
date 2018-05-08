#!/bin/bash

## This file runs the necessary steps to prepare the simulation
## of a given patient(s), entered as input parameters

set -e

# Pull patient from database
# ./_setMCAuto.sh "$@"
# # Once a patient is exported from ERISOne, organize the files and set the weights
# ./_organizeOpt4DFiles.sh "$@"
./_adjustWeights.sh "$@"
# Set the other fraction data and files
# ./_setFractionsFiles.sh "$@"
# ./_setGPUFiles.sh "$@"
