#!/bin/tcsh

# setup environment
setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.0.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.0.0/x directory

set ver = "v2.0.0"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/bc.json  $ver/BC.$ver.json
  ln -s ${PWD}/geosparticles/oc.json  $ver/OC.$ver.json
  ln -s ${PWD}/geosparticles/ni.json  $ver/NI.$ver.json
  ln -s ${PWD}/geosparticles/su.json  $ver/SU.$ver.json
  ln -s ${PWD}/geosparticles/ss.json  $ver/SS.$ver.json
  ln -s ${PWD}/geosparticles/brc.json $ver/BR.$ver.json
  ln -s ${PWD}/geosparticles/du-grasp_spheroid-lognormal.json $ver/DU.$ver.json

# Run the cases

# All wavelengths
  foreach XX (DU BC OC SU BR NI SS)
   ./runoptics.py --name $ver/$XX.$ver.json --dest=$ver> $ver/optics_$XX.$ver.txt &
  end
  wait

# Add phase matrices
  foreach XX (DU BC OC SU BR NI SS)
   ./rungsf.py --filename $ver/optics_$XX.$ver.nomom.nc4 --dest=$ver > $ver/optics_$XX.$ver.gsf.txt &
  end
  wait

# Bands
  foreach XX (DU BC OC SU BR NI SS)
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach XX (DU BC OC SU BR NI SS)
    ./plotoptics.py --name ./AerosolOptics/$ver/x/optics_$XX.$ver.nc4
  end
