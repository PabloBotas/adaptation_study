# Adaptive project scripts and analysis

Scripts and reports of the adaptive proton therapy project

## To add a new patient from Jihun's data
There is more information in the README.md at `patients/reports/README.md`.
Jihun has improved the data formats so this is what it is currently necessary:

- `convertCBCTtofloat.sh`: the code to transform the CBCTs to gPMC format is buggy, but works with floats :)
- `maskSkin.sh`: CAREFUL!! this script will destroy the original CBCTs, use with caution.
- `re-createDicoms.sh`: The dicom directory should be emptied beforehand.

## To add a new patient plan:
- Create directory with patient number and ID in `plans/Opt4D`
- `./setMCAuto.sh` with the directory in the main for loop (edit file or pass as parameter)
- Copy the patient directory (after having run MCAuto) to the cluster and create the Dijs per beam
- Link the Dijs in the Opt4D directory and create the structure.vv files in the local computer
- From ERISOne, send files (there is a script for that) to local computer
- Organize Opt4D files with `_organizeOpt4DFiles.sh` and `_adjustWeights.sh`
- `./setFractionsFiles.sh` with the directory in the main for loop
- `./setGPUFiles.sh`

## To run a plan
The script `runPatient.sh` should run all the patients passed as parameters.
