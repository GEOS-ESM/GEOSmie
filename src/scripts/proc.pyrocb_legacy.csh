#!/bin/tcsh

# setup environment
setenv SRC_DIR /home/colarco/sandbox/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce optical table used in ACMAP23_ANY pyrocb project
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/pyrocb/x directory

set ver = "pyrocb"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/experimental/brc_carma_pyrocb_aged_s12.json \
        $ver/BR.carma_pyrocb_aged_s12.json
  ln -s ${PWD}/geosparticles/experimental/oc_carma_pyrocb_aged_s12.json \
        $ver/OC.carma_pyrocb_aged_s12.json
  ln -s ${PWD}/geosparticles/experimental/brc_carma_pyrocb_aged_s12_bcgf_low3.json \
        $ver/BR.carma_pyrocb_aged_s12_bcgf_low3.json


# Run the cases
  ./runoptics.py -c --name $ver/BR.carma_pyrocb_aged_s12.json \
                 --dest=$ver> $ver/optics_BR.carma_pyrocb_aged_s12.txt &
  ./runoptics.py -c --name $ver/OC.carma_pyrocb_aged_s12.json \
                 --dest=$ver> $ver/optics_OC.carma_pyrocb_aged_s12.txt &
  ./runoptics.py -c --name $ver//BR.carma_pyrocb_aged_s12_bcgf_low3.json \
                 --dest=$ver> $ver/optics_BR.carma_pyrocb_aged_s12_bcgf_low3.txt &
  wait

  ./rungsf.py --filename $ver/optics_BR.carma_pyrocb_aged_s12.nomom.legacy.nc4 --dest=$ver
  ./rungsf.py --filename $ver/optics_OC.carma_pyrocb_aged_s12.nomom.legacy.nc4 --dest=$ver
  ./rungsf.py --filename $ver/optics_BR.carma_pyrocb_aged_s12_bcgf_low3.nomom.legacy.nc4 --dest=$ver


  ./runbands.py --filename $ver/optics_BR.carma_pyrocb_aged_s12.legacy.nc4 --dest=$ver
  ./runbands.py --filename $ver/optics_OC.carma_pyrocb_aged_s12.legacy.nc4 --dest=$ver
  ./runbands.py --filename $ver/optics_BR.carma_pyrocb_aged_s12_bcgf_low3.legacy.nc4 --dest=$ver

# Plot
  mkdir -p plots
  ./plotoptics_legacy.py --name $ver/optics_BR.carma_pyrocb_aged_s12.legacy.nc4
  ./plotoptics_legacy.py --name $ver/optics_OC.carma_pyrocb_aged_s12.legacy.nc4
  ./plotoptics_legacy.py --name $ver/optics_BR.carma_pyrocb_aged_s12_bcgf_low3.legacy.nc4


# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

