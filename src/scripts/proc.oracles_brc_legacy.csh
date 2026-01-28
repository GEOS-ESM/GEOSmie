#!/bin/tcsh

# setup environment
setenv SRC_DIR /home/colarco/sandbox/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce optical table used in ACMAP23_ANY pyrocb project
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/pyrocb/x directory

set ver = "oracles_brc"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/experimental/brc_oracles_bcgf_low3.json \
        $ver/BR.oracles_brc.json

# Run the cases
# BR
  ./runoptics.py -c --name $ver/BR.oracles_brc.json \
                 --dest=$ver> $ver/optics_BR.oracles_brc.txt 
  ./rungsf.py --filename $ver/optics_BR.oracles_brc.nomom.legacy.nc4 --dest=$ver
  ./runbands.py --filename $ver/optics_BR.oracles_brc.legacy.nc4 --dest=$ver

# Plot
  mkdir -p plots
  ./plotoptics_legacy.py --name $ver/optics_BR.oracles_brc.legacy.nc4


# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

