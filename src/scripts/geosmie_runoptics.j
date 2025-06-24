#!/bin/csh -f

#######################################################################
#                     Batch Parameters for Run Job
#######################################################################
 
#SBATCH -J geosmie
#SBATCH --nodes=1
#SBATCH --constraint=mil
#SBATCH --time=01:00:00
#SBATCH -A @GROUPID
#SBATCH -o output_geosaqcgan-%j.log
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --qos=debug
#######################################################################
#  Run GEOS implementation of NASA-AQcGAN to test running a forecast
#######################################################################
#           Architecture Specific Environment Variables
#######################################################################

setenv SRC_DIR @SRCDIR
setenv PYTHONPATH ${SRC_DIR}/install/lib/Python

source $SRC_DIR/env@/g5_modules

#######################################################################
#          Run the GEOSMIE code
# To calculate single-scattering properties at individual wavelengths
# Where file.json defines all of the parameters of the calculation. 
# ".json" can be omitted from the runoptics.py command as a convenience. 
# See JSON examples under geosparticles/ for current GEOS aerosol particles.

#######################################################################

./runoptics.py --name path/to/file.json


