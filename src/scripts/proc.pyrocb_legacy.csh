#!/bin/tcsh

# setup environment
setenv SRC_DIR /gpfsm/dnb06/projects/p233/pcolarco/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.0.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.0.0/x directory

set ver = "pyrocb"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/experimental/brc_carma_pyrocb_aged_s12.json \
        $ver/BR.carma_pyrocb_aged_s12.json

# Run the cases
# BR
  ./runoptics.py -c --name $ver/BR.carma_pyrocb_aged_s12.json \
                 --dest=$ver> $ver/optics_BR.carma_pyrocb_aged_s12.txt 
  ./rungsf.py --filename $ver/optics_BR.carma_pyrocb_aged_s12.nomom.legacy.nc4 --dest=$ver
  ./runbands.py --filename $ver/optics_BR.carma_pyrocb_aged_s12.legacy.nc4 --dest=$ver

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

