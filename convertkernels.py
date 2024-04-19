import os
import numpy as np
import netCDF4 as nc
import json
import glob

def fun(ifn, dest):
  with open(ifn) as fp:
    data = json.load(fp)

  path = data['path']
  fnpre = data['fnpre']
  ratios = data['ratios'] # list of kernel shape id's to use (axis ratios for spheroids)
  contlen = data['contlen'] # length of content (define?)
  mrlen = data['mrlen'] # number of real refractive indices
  milen = data['milen'] # number of imaginary refractive indices
  scathdrlen = data['scathdrlen'] # number of lines in scattering element file headers
  scatcontlen = data['scatcontlen'] # number of content (define?) lines in scattering element files
  elems = data['elems'] # scattering matrix element names
  numang = data['numang'] # number of angle points
  scatelemlen = data['scatelemlen'] # number of lines in scattering matrix element chunks
  kernelname = data['name']

  pfx = path
  sizes, x = getSizes(pfx)

  ratnums = []
  for r in ratios:
    rr = float(r) / 100.
    ratnums.append(rr)
  
  mr = [0] * mrlen
  mi = [0] * milen
  angs = range(numang) # TODO: Does this assume Δθ=1 always? If so, why not numang=180° always?

  allext = []
  allabsorb = []
  allscadata = []

  for rati in range(len(ratios)):
    ratio = ratios[rati]
    print("ri %d of %d, ratio %s"%(rati+1, len(ratios), ratios[rati]))
    exti, absorb, rilist = read00(pfx, ratio, fnpre, contlen)

    # get imag refractive indices
    for i in range(len(mi)):
      mi[i] = rilist[i][1]

    # get real refractive indices
    for i in range(len(mr)):
      mr[i] = rilist[i*len(mi)][0]

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
          # multiplier of 1000 comes from GRASP conventions
          refraind = mri * len(mi) + mii
          thisext = allext[rati][refraind][xi] * 1000. 
          thisabs = allabsorb[rati][refraind][xi] * 1000.

          ext[rati,mri,mii,xi] = thisext

          abso[rati,mri,mii,xi] = thisabs
          for eli in range(len(elems)):
            thissca = allscadata[rati][eli][refraind][xi]
            scama[rati,mri,mii,xi,eli,:] = np.array(thissca) * 1000
          

  sca = ext - abso

  # save everything in netCDF
  print('Opening NetCDF and saving stuff')
  # TODO: Are size parameters saved anywhere?
    # Set output filename and open file
  fn = os.path.join(dest, f"kernel-{kernelname}.nc")
  
  with nc.Dataset(fn, 'w', format='NETCDF4') as ncdf: # open/create netCDF4 data file

    # Create dimensions
    ncdf.createDimension('ratio', len(ratios))
    ncdf.createDimension('mr', len(mr))
    ncdf.createDimension('mi', len(mi))
    ncdf.createDimension('x', len(x))
    ncdf.createDimension('angle', len(angs))
    ncdf.createDimension('scattering_element', len(elems))

    # Create variables
    usezlib = False
    ncdf.createVariable('ratio', 'f8', ('ratio'), zlib=usezlib)
    ncdf.createVariable('mr', 'f8', ('mr'), zlib=usezlib)
    ncdf.createVariable('mi', 'f8', ('mi'), zlib=usezlib)
    ncdf.createVariable('x', 'f8', ('x'), zlib=usezlib)
    ncdf.createVariable('angle', 'f8', ('angle'), zlib=usezlib)

    ncdf.variables['x'][:] = x
    ncdf.variables['ratio'][:] = ratnums

    scalarelems = ('ratio', 'mr', 'mi', 'x')
    scatelems = ('ratio', 'mr', 'mi', 'x', 'scattering_element', 'angle')

    ncdf.createVariable('ext', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('abs', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('sca', 'f8', scalarelems, zlib=usezlib)

    ncdf.createVariable('qext', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('qabs', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('qsca', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('cext', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('cabs', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('csca', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('qb', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('lidar_ratio', 'f8', scalarelems, zlib=usezlib)
    ncdf.createVariable('g', 'f8', scalarelems, zlib=usezlib)

    ncdf.createVariable('scama', 'f8', scatelems, zlib=usezlib)

    ncdf.variables['mr'][:] = mr
    ncdf.variables['mi'][:] = mi
    ncdf.variables['ext'][:,:,:,:] = ext
    ncdf.variables['abs'][:,:,:,:] = abso
    ncdf.variables['scama'][:,:,:,:,:,:] = scama

    # add extra variables
    print('calculate extra variables')

    ncdf.variables['sca'][:,:,:,:] = sca
    volconv = 4./3. * np.array(sizes)

    # grasp kernels need to be divided by a factor of log(x[n+1]/x[n]), i.e. the log of bin size ratios
    # for the size bins used here this ends up being a multiplication by 3.68176736
    graspfactor = 3.68176736 # TODO: Calculate this with above formula from sizes variable.
    ncdf.variables['qsca'][:,:,:,:] = ncdf.variables['sca'][:,:,:,:] * volconv * graspfactor
    ncdf.variables['qext'][:,:,:,:] = ncdf.variables['ext'][:,:,:,:] * volconv * graspfactor

    # calculate cross-sections by multiplying efficiencies by geom cross-section of equivalent spheres
    areaconv = np.pi * np.array(sizes) ** 2
    ncdf.variables['csca'][:,:,:,:] = ncdf.variables['qsca'][:,:,:,:] * areaconv
    ncdf.variables['cext'][:,:,:,:] = ncdf.variables['qext'][:,:,:,:] * areaconv

    f11 = ncdf.variables['scama'][:,:,:,:,0,:]
    f11back = ncdf.variables['scama'][:,:,:,:,0,-1]
    qsca = ncdf.variables['qsca'][:,:,:,:] 
    qext = ncdf.variables['qext'][:,:,:,:] 

    # scama is normalized so we un-normalize it by multiplying by scattering cross-section? <– Below would suggest it is absolute, not normalized...?
    pBck = f11back / sca # get p11 from f11 by dividing by sca
    qBck = pBck * ncdf.variables['qsca'][:,:,:,:] # get backscattering efficiency by multiplying p11 by qsca

    ncdf.variables['qb'][:,:,:,:] = qBck

    angs = np.radians(angs)
    normfactor = np.sum(f11 * np.sin(angs), axis=-1)
    f11n = np.array([f11[:,:,:,:,ai] / normfactor for ai in range(len(angs))])
    g_pre = np.cos(angs) * np.sin(angs) * f11n.T
    g = np.sum(g_pre, axis=-1).T # Note: the omissions of dθ here and in normfactor cancel out
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
  return sizes, x

def readScatEle(pfx, ratio, ele, fnpre, hdrlen, contlen, scaelemlen):
  elename = ele
  elems = readEle(pfx, ratio, hdrlen, contlen, elename, fnpre)
  scama = getScaMaElem(elems, scaelemlen)
  return scama

def read00(pfx, ratio, fnpre, contlen):
  hdrlen = 5 # number of header lines
  #contlen = 16 # how long each "element" is
  elename = '00'
  elems = readEle(pfx, ratio, hdrlen, contlen, elename, fnpre)
  if contlen == 6:
    onelen = 2
  else:
    onelen = 7
  exti, absorb, rilist = getEA(elems, onelen)
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

def getEA(elems, onelen):
  # return separate list of all extinction and all absorb
  exti = []
  absorb = []
  eleheader = 2 # 2 lines for each 
  #foo1 = 7
  foo1 = onelen
  rilist = []
  for ele in elems:
    hdr = ele[:eleheader]
    riline = hdr[1].split()
    mr = riline[1]
    mi = riline[2]
    rilist.append((mr, mi))

    cont = ele[eleheader:]
    thisexti = cont[:foo1][1:]
    thisabs = cont[foo1:][1:]

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
