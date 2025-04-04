#!/usr/bin/env python3.11

import shutil
import os
from optparse import OptionParser
import convertncdf

if __name__ == "__main__":
  parser = OptionParser(usage="Usage: %prog",
                          version='0.0.1' )

  parser.add_option("--filename", dest="filename", default="", 
              help="Optical table file to use (default=%s)"\
                          %(""))

  parser.add_option("--dest", dest="dest", default=".",
              help="Output directory (default=%s)"\
                          %("."))

  parser.add_option("--mode", dest="mode", default="pygeos",
              help="Input file format (default=%s)"\
                          %("pygeos"))

  parser.add_option("--rhop", dest="rhop", default=1000.0,
              help="Particle density (not needed/used for modes pygeos, legendre) (default=%s)"\
                          %(1000.0))

  (options, args) = parser.parse_args()

  if not os.path.exists(options.filename):
    parser.error("Input file path (--filename) does not exist")

  convertncdf.convertFile(options.filename, options.dest, options.mode, options.rhop)
  print('Done!')

