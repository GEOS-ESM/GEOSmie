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

  if not os.path.exists(options.filename):
    parser.error("Input file path (--filename) does not exist")

  if options.bandmode.upper() not in acceptedBandmodes:
    parser.error("Band type must be one of: %s"%(acceptedBandmodes))

  bandaverage.processFileForBandMode(options.filename, options.partname, options.dest,\
  options.bandmode.upper(), options.useSolar, options.noIR)

  print('Done!')

