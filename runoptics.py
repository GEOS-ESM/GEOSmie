#!/usr/bin/env python3

"""
 Utility to create GEOS/GOCART-compatible aerosol optics look-up tables
 Top-level caller for optics table generation code
 Usage example:
  runoptics.py --name bc.json
 will produce an output file named "integ-bc-raw.nc"

  runoptics -h
 will print a usage message

"""

import shutil
import os
import dointegration
import hydrophobic
from optparse import OptionParser
import particleparams as pp

if __name__ == "__main__":
  parser = OptionParser(usage="Usage: %prog",
                          version='0.0.1' )

  acceptedDatatypes = ['json']

  parser.add_option("--name", dest="name", default="", 
              help="Particle file to use (default=%s)"\
                          %(""))

  parser.add_option("--namelist", dest="namelist", default="", 
              help="File with list of particle files (to be passed to --file) to run iteratively. If used, overrides --file (default=%s)"\
                          %(""))
  """
  parser.add_option("--shape", dest="shape", default="mie", 
              help="Particle shape to use %s (default=%s)"\
                          %(acceptedShapetypes, "mie"))
  """

  parser.add_option("--datatype", dest="datatype", default="json", 
              help="Particle data type to use %s (default=%s)"\
                          %(acceptedDatatypes, "json"))

  parser.add_option("--dest", dest="dest", default=".", 
              help="Output directory to use (default=%s)"\
                          %("."))

  parser.add_option("-c","--classic",
                    action="store_true", dest="classic", default=False,
                    help="write output filename is legacy dimensioning")

  # TODO! add flag for running bands either in combination or separately
  # so like mode = [scatonly, bandsonly, both]

  (options, args) = parser.parse_args()

  if options.datatype not in acceptedDatatypes:
    parser.error("data type must be one of: %s"%(acceptedDatatypes))

  if not os.path.exists(options.dest):
    parser.error("Output directory (--dest) does not exist")

  if options.name == "" and options.namelist == "":
    parser.error("non-empty particle name or namelist required (use --name particlename or --namelist namelist)")

  namelist = []
  if options.namelist:
    if not os.path.exists(options.namelist):
      parser.error("Namelist %s does not exist"%options.namelist)

    with open(options.namelist) as fp:
      lines = fp.readlines()

    for line in lines:
      namelist.append(line.strip())
  else:
    namelist = [options.name]

  # Loop over namelisted files
  for fni, fn in enumerate(namelist):
    print("Starting particle file %s, %d of %d"%(fn, fni+1, len(namelist)))
    if not os.path.exists(fn):
      # try with added json in case it was emitted
      fn1 = "%s.json"%fn

      if not os.path.exists(fn1):
        parser.error("File %s doesn't exist"%fn)
      else:
        fn = fn1

    # Given configuration file fn, return particle properties defined
    # Invoked here only to know if special handling for philic/phobic
    # binned files is needed (see later)
    params = pp.getParticleParams(fn, options.datatype)

    # remove path (only use filename) and remove json suffix
    particlename = fn.split('/')[-1].replace(".json", "")

    # Perform the requested calculation/integrations
    dointegration.fun(fn, options.datatype, options.dest, options.classic)

    # Special handling based on key "hydrophobic" in configuration
    # This handling implies a full (hydrophilic) calculation at all
    # RH is done for for one particle bin and the RH=0 (hydrophobic)
    # values are propagated as an additional bin valid at all RH 
    # that is prepended to the results of the initial hydrophilic 
    # calculation
    if options.classic:
      opfn = "integ-%s-raw.legacy.nc"%particlename
    else:
      opfn = "integ-%s-raw.nc"%particlename
    if "hydrophobic" in params and params["hydrophobic"]:
      # rename non-HP file 
      fn2 = "%s.nohp"%opfn
      shutil.move(os.path.join(options.dest, opfn), os.path.join(options.dest, fn2))
      # run the hydrophobic code
      print("Starting hydrophobic bin handling")
      hydrophobic.doConversion(fn2, opfn, options.dest, options.classic)
      # remove the temp file (fn2)
      os.remove(os.path.join(options.dest, fn2))

    print("Done, output file: %s"%opfn)


