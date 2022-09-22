import shutil
import os
from optparse import OptionParser
import bandaverage

if __name__ == "__main__":
  parser = OptionParser(usage="Usage: %prog",
                          version='0.0.1' )

  acceptedBandmodes = ['GEOS5', 'RRTMG', 'RRTMGP', 'PURDUE']

  parser.add_option("--filename", dest="filename", default="", 
              help="Optical table file to use (default=%s)"\
                          %(""))

  parser.add_option("--namelist", dest="namelist", default="", 
              help="File with list of particle optics files (to be passed to --filename) to run iteratively. If used, overrides --filename (default=%s)"\
                          %(""))

  parser.add_option("--partname", dest="partname", default="", 
              help="Particle name to use in the qname variable of the output table (default=%s)"\
                          %(""))

  parser.add_option("--dest", dest="dest", default=".",
              help="Output directory (default=%s)"\
                          %("."))

  parser.add_option("--bandmode", dest="bandmode", default="RRTMG", 
              help="Band averaging type to use %s (default=%s)"\
                          %(acceptedBandmodes, "RRTMG"))

  parser.add_option("--noIR", dest="noIR", default=False, 
              help="Use the noIR option for GEOS5 band type (default=%s)"\
                          %(False))

  parser.add_option("--useSolar", dest="useSolar", default=False, 
              help="Use the useSolar option (default=%s)"\
                          %(False))

  (options, args) = parser.parse_args()

  if options.bandmode.upper() not in acceptedBandmodes:
    parser.error("Band type must be one of: %s"%(acceptedBandmodes))

  namelist = []
  if options.namelist:
    if not os.path.exists(options.namelist):
      parser.error("Namelist %s does not exist"%options.namelist)

    with open(options.namelist) as fp:
      lines = fp.readlines()

    for line in lines:
      namelist.append(line.strip())
  else:
    namelist = [options.filename]

  for fni, fn in enumerate(namelist):
    print("Starting optics file %s, %d of %d"%(fn, fni+1, len(namelist)))

    if not os.path.exists(fn):
      parser.error("Input file path (--filename) does not exist: %s"%fn)
  
    bandaverage.processFileForBandMode(fn, options.partname, options.dest,\
    options.bandmode.upper(), options.useSolar, options.noIR)

  print('Done!')

