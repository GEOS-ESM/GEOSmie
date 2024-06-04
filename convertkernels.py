import os
import numpy as np
import netCDF4 as nc
import json
import glob
import re
from scipy import integrate

def fun(ifn, dest):
  with open(ifn) as fp:
    data = json.load(fp)

  kernelname = data['name'] # name of output file (saved in working directory)
  path = data['path'] # path to the kernel files to be read in
  fnpre = data['fnpre'] # filename preface of kernels to be read in (e.g. "Rkernel1")
  ratios = data['ratios'] # list of kernel shape id's to use (axis ratios for spheroids)
  contlen = data['contlen'] # number of lines per ext+abs block in 00 kernels files (i.e., number of lines between a given "element, ratio" and the next "element, ratio")
  mrlen = data['mrlen'] # number of real refractive indices
  milen = data['milen'] # number of imaginary refractive indices
  scathdrlen = data['scathdrlen'] # number of lines in scattering element file headers
  scatcontlen = data['scatcontlen'] # number of lines for a given [mr,mi] block in phase matrix element kernels files (i.e., number of lines between a given "element, ratio" and the next "element, ratio")
  elems = data['elems'] # 2 digit scattering matrix element indices as strings
  scatelemlen = data['scatelemlen'] # number of lines in scattering matrix element chunks (e.g., one chunk contains 181 angles corresponding to a given (x,n,k))
  xMaxRenorm = data['renormUpperX'] # largest size parameter for which PM elements are renormalized (0 for no renomralization)
  angFwdPeak = data['angFwdPeak'] # angle [degrees] of forward peak in asymmetry parameter calculation truncation correction; 0 -> no correction (e.g., with high angular res. Saito DB) 
  kernelScaleFact = data['kernelScaleFact'] # unity for GRASP v1.1.3 kernels; kernel values multiplied by kernelScaleFact should give units of 1/μm 

  pfx = path
  sizes, x, lamb = getSizes(pfx)

  ratnums = [float(r)/100 for r in ratios]

  allext = []
  allabsorb = []
  allscadata = []

  for rati in range(len(ratios)):
    ratio = ratios[rati]
    print("ri %d of %d, ratio %s"%(rati+1, len(ratios), ratios[rati]))
    exti, absorb, rilist = read00(pfx, ratio, fnpre, contlen)

    # get refractive indices
    mi = [rilist[i][1] for i in range(milen)]
    mr = [rilist[i*milen][0] for i in range(mrlen)]

    angs = readScaAngles(pfx, ratio, fnpre, scathdrlen)
    scadata = []
    for ei, ele in enumerate(elems):
      thisscadata = readScatEle(pfx, ratio, ele, fnpre, scathdrlen, scatcontlen, scatelemlen)
      scadata.append(thisscadata)
    
    allext.append(exti)
    allabsorb.append(absorb)
    allscadata.append(scadata)

  ext = np.zeros([len(ratios), len(mr), len(mi), len(x)])
  abso = np.zeros_like(ext)

  scama = np.zeros([len(ratios), len(mr), len(mi), len(x), len(elems), len(angs)])
  
  for rati in range(len(ratios)):
    print('rati %d of %d'%(rati+1, len(ratios)))
    for mri in range(len(mr)):
      for mii in range(len(mi)):
        for xi in range(len(x)):
          # calculate the 1d index for refractive index array
          refraind = mri * len(mi) + mii
          thisext = allext[rati][refraind][xi] * kernelScaleFact 
          thisabs = allabsorb[rati][refraind][xi] * kernelScaleFact
          ext[rati,mri,mii,xi] = thisext
          abso[rati,mri,mii,xi] = thisabs
          for eli in range(len(elems)):
            thissca = allscadata[rati][eli][refraind][xi]
            scama[rati,mri,mii,xi,eli,:] = np.array(thissca) * kernelScaleFact
     
     
  # clean up units and add extra variables
  print('calculating extra variables')     
  
  # grasp kernels need to be divided by a factor of log(x[n+1]/x[n]), i.e. the log of bin size ratios
  graspfactors = 1/np.log(np.divide(x[1:],x[:-1])) # this is separate from from "kernelScaleFact" above
  assert np.any(np.diff(graspfactors) < 1e-4), 'Size bins in GRASP kernels must be log-spaced.'
  graspfactor = graspfactors.mean() # ~3.681765 for GRASP v1.1.3 standard kernels
  ext = ext * graspfactor
  abso = abso * graspfactor
  sca = ext - abso
  sca[sca==0] = abso[sca==0]/1e7 # sca==0 is nonphysical but ext and abso only stored to 7 decimal places in the kernel text files
  scama =  (scama*graspfactor)/sca[:,:,:,:,None,None]

  # Phase matrix element renormilzation 
  ang_rad = np.radians(angs)
  kern = np.sin(ang_rad)
  if xMaxRenorm>0: # Renormalize below x<xMaxRenorm (GRASP normalization noisy at x<0.1)
    lowxInd = np.asarray(x) < xMaxRenorm
    p11_raw = scama[:,:,:,lowxInd,0,:]
    nf1 = integrate.simpson(p11_raw * kern, x=ang_rad, axis=-1)/2
    # Renormalize all PM elements to preserve their ratios which appear correct
    scama[:,:,:,lowxInd,:,:] = scama[:,:,:,lowxInd,:,:]/nf1[:,:,:,:,None,None]
  
  # Integrate over scattering angle to find asymmetry parameters
  noPeaki = np.asarray(angs) >= angFwdPeak # define region outside of forward peak for which integration is assumed accurate
  if angFwdPeak > 0: # if >20 angles exist within forward peak (e.g., in Saito DB) we assume traditional integration is accurate 
    nf_noPeak = integrate.simpson(scama[:,:,:,:,0,noPeaki] * kern[noPeaki], x=ang_rad[noPeaki], axis=-1)
    nf_peak = 2 - nf_noPeak # the contribution from the peak summed with above integral should give 2 since P11 is in units of sr-1
  else: # we will determine g by integrating over the full range [0°,180°]
    nf_peak = 0 # there is no separate contribution from the forward peak because it will be included in the main integral
  kern_g = np.cos(ang_rad[noPeaki]) * kern[noPeaki] # asymmetry parameter is integral of p11 weighted by cos(θ)sin(θ)
  g_noPeak = integrate.simpson(scama[:,:,:,:,0,noPeaki] * kern_g, x=ang_rad[noPeaki], axis=-1)
  g = (g_noPeak + nf_peak)/2 # cos(θ)sin(θ)->sin(θ) as θ->0° so the contribution to nf from the forward peak approximates its contribution to g
  
  # Calculate optical efficiencies 
  volconv = 4./3. * np.array(sizes)
  qext = ext * volconv
  qsca = sca * volconv

  # Calculate backscattering efficiencies 
  p11back = scama[:,:,:,:,0,-1]
  qBck = p11back * qsca # get backscattering efficiency by multiplying p11 by qsca

  # calculate cross-sections by multiplying efficiencies by geom cross-section of equivalent spheres
  areaconv = np.pi * np.array(sizes) ** 2
  csca = qsca * areaconv
  cext = qext * areaconv
  
  
  # save everything in netCDF
  print('Opening NetCDF and saving stuff')
  # Set output filename and open file
  fn = os.path.join(dest, f"kernel-{kernelname}.nc")
  
  with nc.Dataset(fn, 'w', format='NETCDF4') as ncdf: # open/create netCDF4 data file

    # Create dimensions
    ncdf.createDimension('ratio', len(ratios))
    ncdf.createDimension('mr', len(mr))
    ncdf.createDimension('mi', len(mi))
    ncdf.createDimension('x', len(x))
    ncdf.createDimension('angle', len(angs))
    ncdf.createDimension('scattering_element_index', len(elems))

    # Create 1D variables
    usezlib = False
    ncdf.createVariable('ratio', 'f8', ('ratio'), zlib=usezlib)
    ncdf.variables['ratio'].long_name = 'aspect ratio'
    ncdf.createVariable('mr', 'f8', ('mr'), zlib=usezlib)
    ncdf.variables['mr'].long_name = 'real refractive index'
    ncdf.createVariable('mi', 'f8', ('mi'), zlib=usezlib)
    ncdf.variables['mi'].long_name = 'imaginary refractive index'
    ncdf.createVariable('x', 'f8', ('x'), zlib=usezlib)
    ncdf.variables['x'].long_name = 'size parameter (2π*r/λ)'
    ncdf.createVariable('angle', 'f8', ('angle'), zlib=usezlib)
    ncdf.variables['angle'].long_name = 'scattering angle'
    ncdf.variables['angle'].units = 'degree'
    ncdf.createVariable('scattering_element_index', 'u1', ('scattering_element_index'), zlib=usezlib)

    # Create Multidimensional variables
    scalarelems = ('ratio', 'mr', 'mi', 'x')
    scatelems = ('ratio', 'mr', 'mi', 'x', 'scattering_element_index', 'angle')
    desc = "This is extinction cross section per unit of particle volume.\
       Alternatively, the extinction coefficient for a volume concentration of unity.\
       ext * ρ = βext where ρ is particle density and βext is mass extinction efficiency."
    varDict = {'ext': 'extinction', 'abs': 'absorption', 'sca': 'scattering'}
    for short,full in varDict.items():
      longname = 'volume %s efficiency' % full
      descNow = desc.replace('extinction',full).replace('ext',short)
      nc4VarSetup(ncdf, short, '1/um', scalarelems, longname, desc=descNow)
      longname = '%s efficiency' % full
      nc4VarSetup(ncdf, 'q'+short, 'none', scalarelems, longname)
      longname = '%s crosss section at λ=%5.3fμm' % (full, lamb)
      nc4VarSetup(ncdf, 'c'+short, 'μm^2', scalarelems, longname)
    nc4VarSetup(ncdf, 'qb', '1/sr', scalarelems, 'backscattering efficiency')
    nc4VarSetup(ncdf, 'lidar_ratio', '1/sr', scalarelems, 'lidar ratio')
    nc4VarSetup(ncdf, 'g', 'none', scalarelems, 'asymmetry parameter')
    desc = 'Normalized such that p11(θ)*sin(θ)*dθ intgrated from 0 to π equals 2.'
    nc4VarSetup(ncdf, 'scama', '1/sr', scatelems, 'scattering matrix element values', desc=desc)

    # Write data to variables
    ncdf.variables['x'][:] = x
    ncdf.variables['ratio'][:] = ratnums
    ncdf.variables['mr'][:] = mr
    ncdf.variables['mi'][:] = mi
    ncdf.variables['angle'][:] = angs
    ncdf.variables['scattering_element_index'][:] = [int(el) for el in elems]
    ncdf.variables['ext'][:,:,:,:] = ext
    ncdf.variables['abs'][:,:,:,:] = abso
    ncdf.variables['sca'][:,:,:,:] = sca
    ncdf.variables['scama'][:,:,:,:,:,:] = scama
    ncdf.variables['qsca'][:,:,:,:] = qsca
    ncdf.variables['qext'][:,:,:,:] = qext
    ncdf.variables['qabs'][:,:,:,:] = qext - qsca
    ncdf.variables['csca'][:,:,:,:] = csca
    ncdf.variables['cext'][:,:,:,:] = cext
    ncdf.variables['cabs'][:,:,:,:] = cext - csca
    ncdf.variables['qb'][:,:,:,:] = qBck
    ncdf.variables['g'][:,:,:,:] = g


def getSizes(pfx):
  grid = 'grid1.*'
  fn = os.path.join(pfx, grid)
  fn = glob.glob(fn)[0] # extension varies (grid1.dat or grid1.txt)
  with open(fn) as fp:
    lines = fp.readlines()
    lines = [line.rstrip() for line in lines]
  hdr = lines[0]
  numsizes = int(hdr.split()[0])
  lamb = float(hdr.split()[1])
  sizes = lines[1:numsizes+1]
  sizes = [float(s) for s in sizes]
  x = [s / lamb * 2 * np.pi for s in sizes] # size parameter
  return sizes, x, lamb


def readScatEle(pfx, ratio, ele, fnpre, hdrlen, contlen, scaelemlen):
  elename = ele
  elems = readEle(pfx, ratio, hdrlen, contlen, elename, fnpre)
  scama = getScaMaElem(elems, scaelemlen)
  return scama

def read00(pfx, ratio, fnpre, contlen):
  hdrlen = 5 # number of header lines
  elename = '00'
  elems = readEle(pfx, ratio, hdrlen, contlen, elename, fnpre)
  if contlen == 6:
    onelen = 2
  else:
    onelen = 7
  exti, absorb, rilist = getEA(elems, contlen)
  return exti, absorb, rilist


def readEle(pfx, ratio, hdrlen, contlen, elename, fnpre):
  fn0 = '%s_%s_%s.txt'%(fnpre, ratio, elename)
  fn = os.path.join(pfx, fn0)
 
  with open(fn) as fp:
    lines = fp.readlines()
    lines = [line.rstrip() for line in lines]

  header = lines[:hdrlen]
  txtdata = lines[hdrlen:]

  numelemline = header[-1].split()[:2] # number of refractive indices we use
  numelem1 = int(numelemline[0]) * -int(numelemline[1])

  numelem2 = len(txtdata) / contlen
  elems = []
  for eli in range(numelem1):
    start = eli * contlen
    end = (eli+1) * contlen
    elems.append(txtdata[start:end])
  
  return elems


def readScaAngles(pfx, ratio, fnpre, hdrlen):
  fn0 = '%s_%s_%s.txt'%(fnpre, ratio, '11')
  fn = os.path.join(pfx, fn0)
  with open(fn) as fp:
    header = [next(fp) for _ in range(hdrlen)]
    header = [line.rstrip() for line in header]
    
  mtchPtrn = '[ \t]*([0-9]+)[ ]*number of (?:scattering )?angles'
  scatAngMatch = [re.match(mtchPtrn, line) for line in header]
  assert np.any(scatAngMatch), 'Scattering angles could not be parsed in header of scattering matrix element files.'
  scatAngLn = np.nonzero(scatAngMatch)[0][0] # line number with scattering angle header (e.g., ' 181   number of scattering angles')
  nScatAngs = int(scatAngMatch[scatAngLn].group(1)) # the number of scattering angles
  scatAngsStr = ' '.join(header[scatAngLn+1:]).split()[:nScatAngs] # the scattering angles themselves (as strings)
  scatAngsFloat = [float(sca) for sca in scatAngsStr]
  
  return scatAngsFloat
  
  
def getScaMaElem(elems, onelen):
  # return separate list of all scattering matrix values
  scama = []
  eleheader = 2 # 2 lines for each 
  for ele in elems:
    cont = ele[eleheader:]
    totlen = len(ele)-eleheader
    numitems = totlen/onelen
    elescama = []
    for i in range(int(numitems)):
      startind = i*onelen
      endind = (i+1)*onelen
      thiscont = cont[startind:endind]
      thisscama = getOneScaMa(thiscont)
      elescama.append(thisscama)
    scama.append(elescama)

  return scama
  
  
def getOneScaMa(li):
  ret = []
  for e in li:
    ret += e.split()
  ret = [float(x) for x in ret]
  return ret


def getEA(elems, contlen):
  # return separate list of all extinction and all absorb
  exti = []
  absorb = []
  eleheader = 2 # each element begins with 2 header lines (element, ratio and wavel, rreal, rimag) 
  assert not (contlen-eleheader) % 2, 'Number of data lines in block can not be odd. Is contlen correct?'
  blockSize = int((contlen-eleheader)/2)
  rilist = []
  for ele in elems:
    hdr = ele[:eleheader]
    riline = hdr[1].split()
    mr = riline[1]
    mi = riline[2]
    rilist.append((mr, mi))

    contnt = ele[eleheader:]
    thisexti = contnt[:blockSize][1:] # indexing [1:] removes "EXTINCTION" or "ABSOPRTION" header line
    thisabs = contnt[blockSize:][1:]

    theseexti = getEAOne(thisexti)
    theseabs = getEAOne(thisabs)
    
    exti.append(theseexti)
    absorb.append(theseabs)

  return exti, absorb, rilist


def getEAOne(li):
  ret = []
  for e in li:
    ret += e.split()
  ret = [float(x) for x in ret]
  return ret


def nc4VarSetup(ncdf, nc4Key, units, dims, lngNm=None, varTyp='f8', desc=None, usezlib=False):
  varHnd = ncdf.createVariable(nc4Key, varTyp, dims, zlib=usezlib)
  varHnd.units = units
  if desc: varHnd.description = desc
  if lngNm: varHnd.long_name = lngNm
  
