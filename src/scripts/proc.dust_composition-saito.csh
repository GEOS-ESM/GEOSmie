#!/bin/tcsh

# setup environment
setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.2.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.2.0/x directory

set ver = "log"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

set sites = ("Algeria" "Arizona" "Atacama" "Australia" \
         "Bodele" "Ethiopia" "Gobi" "Kuwait" \
         "Libya" "Mali" "Mauritania" "Morocco" \
         "Namib-1" "Namib-2" "Niger" "Patagonia" \
         "SaudiArabia" "Taklimakan" "Tunisia")

# Link the desired files
  foreach site (`echo $sites`)
  sed 's/SITE/'"$site"'/g'  ${PWD}/geosparticles/experimental/dust_composition/du-lognorm-saito.json > $ver/DU_saito_$site.$ver.json
  end

# Run the cases

# All wavelengths
  foreach site (`echo $sites`)
   ./runoptics.py --name $ver/DU_saito_$site.$ver.json --dest=$ver> $ver/optics_DU_saito_$site.$ver.txt &
  end
  wait

# Add phase matrices
#  foreach XX (DU_grasp_colarco DU_grasp_woodward DU_grasp_balkanski DU_grasp_balkanski27 DU_grasp_opac\
#              DU_mie_colarco DU_mie_woodward DU_mie_balkanski DU_mie_balkanski27 DU_mie_opac)
#   ./rungsf.py --filename $ver/optics_$XX.$ver.nomom.nc4 --dest=$ver > $ver/optics_$XX.$ver.gsf.txt &
#  end
#  wait

# Bands
#  foreach XX (DU_grasp_colarco DU_grasp_woodward DU_grasp_balkanski DU_grasp_balkanski27 DU_grasp_opac\
#              DU_mie_colarco DU_mie_woodward DU_mie_balkanski DU_mie_balkanski27 DU_mie_opac)
#   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
#   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver --bandmode=RRTMGP
#  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach site (`echo $sites`)
   ./plotoptics_bins.py --name ./AerosolOptics/$ver/x/optics_DU_saito_$site.$ver.nomom.nc4
  end
