# full suite assumes geosparticles/du.json has a valid path to dust kernels
# which are not included by default in GEOSmie distribution
#./runoptics.py --namelist namelists/geos.txt
#./runbands.py --namelist namelists/geos_bands.txt

# following namelists exclude dust and can be run directly after installing GEOSmie
./runoptics.py --dest x --namelist namelists/geos_nodust.txt
./runbands.py --dest x --namelist namelists/geos_nodust_bands.txt
