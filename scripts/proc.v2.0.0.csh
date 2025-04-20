#!/bin/tcsh

# Script to reproduce v2.0.0 optics tables
# JSON files live in geosparticles
# Output will be placed in AerosolOptics/v2.0.0/x directory

set ver = "v2.0.0"

mkdir $ver
mkdir -p ../AerosolOptics/$ver/x

# Link the desired files
  ln -s ../../geosparticles/bc.json  $ver/BC.$ver.json
  ln -s ../../geosparticles/oc.json  $ver/OC.$ver.json
  ln -s ../../geosparticles/ni.json  $ver/NI.$ver.json
  ln -s ../../geosparticles/su.json  $ver/SU.$ver.json
  ln -s ../../geosparticles/ss.json  $ver/SS.$ver.json
  ln -s ../../geosparticles/brc.json $ver/BR.$ver.json
  ln -s ../../geosparticles/du-grasp_spheroid-lognormal.json $ver/DU.$ver.json

# Copy the relevant drivers here
  setenv PYTHONPATH /home/colarco/GEOSmie:/home/colarco/GEOSmie/gsf:$PYTHONPATH
  ln -s ../data data
  ln -s ../kernels kernels
  \cp -f ../runbands.py .
  \cp -f ../runoptics.py .
  \cp -f ../gsf/rungsf.py .
  \cp -f ../gsf/a.out .

# Run the cases

# All wavelengths
  foreach XX (DU BC OC SU BR NI SS)
   ./runoptics.py --name $ver/$XX.$ver.json --dest=$ver> $ver/optics_$XX.$ver.txt &
  end
  wait

# Add phase matrices
  foreach XX (DU BC OC SU BR NI SS)
   ./rungsf.py --filename $ver/optics_$XX.$ver.nomom.nc4 --dest=$ver
  end

# Bands
  foreach XX (DU BC OC SU BR NI SS)
   ./runbands.py --filename $ver/optics_$XX.$ver.nc4 --dest=$ver
  end

# Move files
  \mv -f $ver/*nc4 ../AerosolOptics/$ver/x


# Clean up
#  rm -f runbands.py runoptics.py rungsf.py a.out

