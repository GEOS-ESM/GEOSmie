#!/bin/tcsh

# setup environment
setenv SRC_DIR /gpfsm/dnb06/projects/p233/pcolarco/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce Hunga volcanic sulfate optical table used in APARC experiments
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/hunga/x directory

set ver = "hunga"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/experimental/su_ht_reff04_sig16_single.json
        $ver/SU.hunga_reff04_sig16_single.json

# Run the cases
# SU
  ./runoptics.py -c --name $ver/SU.hunga_reff04_sig16_single.json \
                 --dest=$ver> $ver/optics_SU.hunga_reff04_sig16_single.txt
  ./rungsf.py --filename $ver/optics_SU.hunga_reff04_sig16_single.nomom.legacy.nc4 --dest=$ver
  ./runbands.py --filename $ver/optics_SU.hunga_reff04_sig16_single.legacy.nc4 --dest=$ver

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

