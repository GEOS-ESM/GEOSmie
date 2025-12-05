#!/bin/tcsh

# setup environment
setenv SRC_DIR /home/colarco/sandbox/GEOSmie
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.1.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.1.0/x directory

set ver = "v2.1.0"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
# Create the default CARMA json files
  ./carma_utils.py --species SU --rlow 2.0e-10 --rup 7.84672102e-6 --versionid $ver
  \mv -f carmaSU.$ver.json carmaSUpure.$ver.json
  ./carma_utils.py --species SU --versionid $ver
  ./carma_utils.py --species DU --rhop 2650.0 --versionid $ver
  ./carma_utils.py --species SS --rhop 2200.0 --versionid $ver
  ./carma_utils.py --species OC --rhop 1800.0 --versionid $ver
  ./carma_utils.py --species BC --rhop 1800.0 --versionid $ver

# Run the cases

# All wavelengths
  foreach XX (DU BC OC SU SS SUpure)
   ./runoptics.py -c --name carma$XX.$ver.json --dest=$ver> $ver/optics_carma$XX.$ver.txt &
  end
  wait

# Add phase matrices
  foreach XX (DU BC OC SU SS SUpure)
   ./rungsf.py --filename $ver/optics_carma$XX.$ver.nomom.legacy.nc4 --dest=$ver> $ver/optics_carma$XX.$ver.gsf.txt &
  end
  wait

# Bands
  foreach XX (DU BC OC SU SS SUpure)
   ./runbands.py --filename $ver/optics_carma$XX.$ver.legacy.nc4 --dest=$ver
  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach XX (DU BC OC SU SS SUpure)
    ./plotoptics_bins_legacy.py --name ./AerosolOptics/$ver/x/optics_carma$XX.$ver.legacy.nc4
  end
  \mv -f *png plots

