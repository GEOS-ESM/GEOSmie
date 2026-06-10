#!/bin/tcsh

# setup environment
setenv SRC_DIR /gpfsm/dnb06/projects/p233/pcolarco/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to calculate sulfate optical properties
# Case: based on particle property assumptions in Spurr et al. 2026
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/sulfate/x directory

set ver = "sulfate"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/experimental/su_spurr_2026.json \
        $ver/SU.spurr_2026.json

# Run the cases
# SU
  ./runoptics.py -c --name $ver/SU.spurr_2026.json \
                 --dest=$ver> $ver/optics_SU.spurr_2026.txt
  ./rungsf.py --filename $ver/optics_SU.spurr_2026.nomom.legacy.nc4 --dest=$ver
  ./runbands.py --filename $ver/optics_SU.spurr_2026.legacy.nc4 --dest=$ver

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  ./plotoptics.py --name ./AerosolOptics/$ver/x/optics_SU.spurr_2026.nc4
