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

  # TODO! add flag for running bands either in combination or separately
  # so like mode = [scatonly, bandsonly, both]

  (options, args) = parser.parse_args()

  if options.datatype not in acceptedDatatypes:
    parser.error("data type must be one of: %s"%(acceptedDatatypes))

  if not os.path.exists(options.dest):
    parser.error("Output directory (--dest) does not exist")

  if options.name == "":
    parser.error("non-empty particle name required (use --name particlename)")

  fn = "%s"%options.name

  if not os.path.exists(fn):
    # try with added json in case it was emitted
    fn1 = "%s.json"%fn

    if not os.path.exists(fn1):
      parser.error("File %s doesn't exist"%fn)
    else:
      fn = fn1

  params = pp.getParticleParams(fn, options.datatype)

  # remove path (only use filename) and remove json suffix
  particlename = fn.split('/')[-1].replace(".json", "")

  dointegration.fun(fn, options.datatype, options.dest)


  opfn = "integ-%s-raw.nc"%particlename
  if "hydrophobic" in params and params["hydrophobic"]:
    # rename non-HP file 
    fn2 = "%s.nohp"%opfn
    shutil.move(os.path.join(options.dest, opfn), os.path.join(options.dest, fn2))
    # run the hydrophobic code
    print("Starting hydrophobic bin handling")
    hydrophobic.doConversion(fn2, opfn, options.dest)
    # remove the temp file (fn2)
    os.remove(os.path.join(options.dest, fn2))

  print("Done, output file: %s"%opfn)


