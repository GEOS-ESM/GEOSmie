#!/bin/tcsh

# setup environment
setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.2.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.2.0/x directory

set ver = "dndr"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-colarco.json  $ver/DU_colarco.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-woodward.json  $ver/DU_woodward.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-balkanski.json  $ver/DU_balkanski.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-balkanski27.json  $ver/DU_balkanski27.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-opac.json  $ver/DU_opac.$ver.json

# Run the cases

# All wavelengths
  foreach XX (DU_colarco DU_woodward DU_balkanski DU_balkanski27 DU_opac)
   ./runoptics.py --name $ver/$XX.$ver.json --dest=$ver> $ver/optics_$XX.$ver.txt &
  end
  wait

# Add phase matrices
  foreach XX (DU_colarco DU_woodward DU_balkanski DU_balkanski27 DU_opac)
   ./rungsf.py --filename $ver/optics_$XX.$ver.nomom.nc4 --dest=$ver > $ver/optics_$XX.$ver.gsf.txt &
  end
  wait

# Bands
  foreach XX (DU_colarco DU_woodward DU_balkanski DU_balkanski27 DU_opac)
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver --bandmode=RRTMGP
  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach XX (DU_colarco DU_woodward DU_balkanski DU_balkanski27 DU_opac)
    ./plotoptics.py --name ./AerosolOptics/$ver/x/optics_$XX.$ver.nc4
  end
