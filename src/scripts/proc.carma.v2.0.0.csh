#!/bin/tcsh

# setup environment
setenv SRC_DIR /gpfsm/dnb34/pcolarco/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.0.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.0.0/x directory

set ver = "v2.0.0"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
# Create the default CARMA json files
  foreach XX (DU BC OC SU SS)
   ./carma_utils.py --species $XX
  end
  
# Run the cases

# All wavelengths
  foreach XX (DU BC OC SU SS)
   ./runoptics.py -c --name carma$XX.json --dest=$ver> $ver/optics_carma$XX.$ver.txt &
  end
  wait

# Add phase matrices
  foreach XX (DU BC OC SU SS)
   ./rungsf.py --filename $ver/optics_carma$XX.$ver.nomom.nc4 --dest=$ver
  end

# Bands
  foreach XX (DU BC OC SU SS)
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach XX (DU BC OC SU SS)
    ./plotoptics_bins_legacy.py --name ./AerosolOptics/$ver/x/optics_carma$XX.$ver.nc4
  end
