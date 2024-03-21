# full suite assumes geosparticles/du.json has a valid path to dust kernels
# which are not included by default in GEOSmie distribution
#python3.8 runoptics.py --namelist namelists/geos.txt
#python3.8 runbands.py --namelist namelists/geos_bands.txt

# following namelists exclude dust and can be run directly after installing GEOSmie
python3.8 runoptics.py --namelist namelists/geos_nodust.txt
python3.8 runbands.py --namelist namelists/geos_nodust_bands.txt
