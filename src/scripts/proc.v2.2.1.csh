#!/bin/tcsh
#SBATCH --time=12:00:00
##SBATCH --nodes=1 --ntasks-per-node=125
#SBATCH --ntasks=12
#SBATCH --job-name=GEOSmie_v2.2.1
#SBATCH --account=s2191

# setup environment
setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

# Script to reproduce v2.2.1 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.2.1/x directory

set ver = "v2.2.1"

mkdir $ver
mkdir -p ./AerosolOptics/$ver/x

# Link the desired files
  ln -s ${PWD}/geosparticles/bc.json  $ver/BC.$ver.json
  ln -s ${PWD}/geosparticles/oc.json  $ver/OC.$ver.json
  ln -s ${PWD}/geosparticles/ni.json  $ver/NI.$ver.json
  ln -s ${PWD}/geosparticles/su.json  $ver/SU.$ver.json
  ln -s ${PWD}/geosparticles/ss.json  $ver/SS.$ver.json
  ln -s ${PWD}/geosparticles/experimental/brc_oracles_bcgf_low3.json $ver/BR.$ver.json
  ln -s ${PWD}/geosparticles/du-grasp_spheroid.json $ver/DU.$ver.json

# Run the cases

# All wavelengths
  set json = `\ls -1 $ver/*.$ver.json`
  prund.pl -H `hostname` -d `echo $json` &
  mpirun -np 12 prund.pl -H `hostname` ./runoptics.py --name %s --dest=$ver> $ver/optics_$ver.txt
  wait

# Add phase matrices
  set infile = `\ls -1 $ver/*.$ver.nomom.nc4`
  prund.pl -H `hostname` -d `echo $infile` &
  mpirun -np 12 ./rungsf.py --filename %s --dest=$ver > $ver/optics_$ver.gsf.txt
  wait

# Bands
  foreach XX (DU BC OC SU BR NI SS)
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver --bandmode=RRTMGP
  end

# Move files
  \mv -f $ver/*nc4 ./AerosolOptics/$ver/x

# Make plots
  mkdir -p plots
  foreach XX (DU BC OC SU BR NI SS)
    ./plotoptics.py --name ./AerosolOptics/$ver/x/optics_$XX.$ver.nc4
  end
