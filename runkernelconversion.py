#!/usr/bin/env python3.11

"""
 Utility to create GEOSmie compatible netcdf kernel 
 kernel files from input GRASP-like text files.
 Usage example:
  runkernelconversion.py --filename data/kernelconversion/grasp.json
 will produce an output file named "kernel-grasp-v1.1.3-integrated_V4.nc"

"""
import shutil
import os
from optparse import OptionParser
import convertkernels

if __name__ == "__main__":
  parser = OptionParser(usage="Usage: %prog",
                          version='0.0.1' )

  parser.add_option("--filename", dest="filename", default="", 
              help="Kernel conversion JSON file path (default=%s)"\
                          %(""))

  parser.add_option("--dest", dest="dest", default=".",
              help="Output directory (default=%s)"\
                          %("."))

  (options, args) = parser.parse_args()

  if not os.path.exists(options.filename):
    parser.error("Input file path (--filename) does not exist")

  convertkernels.fun(options.filename, options.dest)
  print('Done!')

