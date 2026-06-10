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
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-colarco.json  $ver/DU_grasp_colarco.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-woodward.json  $ver/DU_grasp_woodward.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-balkanski.json  $ver/DU_grasp_balkanski.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-balkanski27.json  $ver/DU_grasp_balkanski27.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-grasp_spheroid-opac.json  $ver/DU_grasp_opac.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-mie-colarco.json  $ver/DU_mie_colarco.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-mie-woodward.json  $ver/DU_mie_woodward.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-mie-balkanski.json  $ver/DU_mie_balkanski.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-mie-balkanski27.json  $ver/DU_mie_balkanski27.$ver.json
  ln -s ${PWD}/geosparticles/experimental/dust/du-mie-opac.json  $ver/DU_mie_opac.$ver.json

# Run the cases

# All wavelengths
  foreach XX (DU_grasp_colarco DU_grasp_woodward DU_grasp_balkanski DU_grasp_balkanski27 DU_grasp_opac\
              DU_mie_colarco DU_mie_woodward DU_mie_balkanski DU_mie_balkanski27 DU_mie_opac)
   ./runoptics.py --name $ver/$XX.$ver.json --dest=$ver> $ver/optics_$XX.$ver.txt &
  end
  wait

# Add phase matrices
  foreach XX (DU_grasp_colarco DU_grasp_woodward DU_grasp_balkanski DU_grasp_balkanski27 DU_grasp_opac\
              DU_mie_colarco DU_mie_woodward DU_mie_balkanski DU_mie_balkanski27 DU_mie_opac)
#   ./rungsf.py --filename $ver/optics_$XX.$ver.nomom.nc4 --dest=$ver > $ver/optics_$XX.$ver.gsf.txt &
  end
  wait

# Bands
  foreach XX (DU_grasp_colarco DU_grasp_woodward DU_grasp_balkanski DU_grasp_balkanski27 DU_grasp_opac\
              DU_mie_colarco DU_mie_woodward DU_mie_balkanski DU_mie_balkanski27 DU_mie_opac)
#   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
#   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver --bandmode=RRTMGP
  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach XX (DU_grasp_colarco DU_grasp_woodward DU_grasp_balkanski DU_grasp_balkanski27 DU_grasp_opac\
              DU_mie_colarco DU_mie_woodward DU_mie_balkanski DU_mie_balkanski27 DU_mie_opac)
   ./plotoptics.py --name ./AerosolOptics/$ver/x/optics_$XX.$ver.nmom.nc4
  end
