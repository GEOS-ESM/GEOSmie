import os
import numpy as np
import netCDF4 as nc

def fun(mode):
  if mode == 'grasp':
    path = os.path.join('/misc', 'opk01', 'spheroid-package-master')
    pfx0 = 'KERNELS_n22k16_181_base' 
    fnpre = 'Rkernel1'
    ratios = ['033', '037', '040', '044', '048', '053', '058',\
    '063', '069', '076', '083', '091', '100', '109', '120', '131',\
    '144', '158', '173', '189', '207', '227', '249', '273', '299']
    opn = 'test15-temp.nc'
    contlen = 16
    mrlen = 22
    milen = 16
    scathdrlen = 32 # number of header lines
    scatcontlen = 1068 # how long each "element" is
    #elems = ['11', '12', '22', '33', '34', '44']
    elems = ['11', '22', '33', '44', '12', '34'] # this is the assumed read in order
    numang = 181
    scatelemlen = 26 # how many lines per one refra
  elif mode == 'saito':
    path = os.path.join('/misc', 'opk01', 'saito-kernels')
    pfx0 = ''
    fnpre = 'kernels'
    ratios = ['299']
    opn = 'saito-test.nc'
    contlen = 6 # no hard line wraps unlike grasp
    mrlen = 8
    milen = 31
    scathdrlen = 25 # saito
    scatcontlen = 275 # saito
    #elems = ['11', '12', '22', '33', '43', '44']
    #elems = ['11', '12', '33', '43', '22', '44']
    elems = ['11', '22', '33', '44', '12', '43'] # this is the assumed read in order
    numang = 4706
    numang = 181
    scatelemlen = 1

  pfx = os.path.join(path, pfx0)
  sizes, x = getSizes(pfx)

  ratnums = []
  for r in ratios:
    rr = float(r) / 100.
    ratnums.append(rr)
  
  # output filename
  fn = os.path.join('../../opk01/convertspheroids/', opn)

  ncdf = nc.Dataset(fn, 'w')
  mr = [0] * mrlen
  mi = [0] * milen
  angs = range(numang)

  # Create dimensions
  ncdf.createDimension('ratio', len(ratios))
  ncdf.createDimension('mr', len(mr))
  ncdf.createDimension('mi', len(mi))
  ncdf.createDimension('x', len(x))
  ncdf.createDimension('angle', len(angs))
  ncdf.createDimension('scattering_element', len(elems))

  usezlib = False

  # Create variables
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

  ncdf.variables['mr'][:] = mr
  ncdf.variables['mi'][:] = mi

  ext = np.zeros([len(ratios), len(mr), len(mi), len(x)])
  abso = np.zeros_like(ext)

  scama = np.zeros([len(ratios), len(mr), len(mi), len(x), len(elems), len(angs)])

  # save everything in netCDF
  print('start netcdf save stuff')
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
  ncdf.variables['ext'][:,:,:,:] = ext
  ncdf.variables['abs'][:,:,:,:] = abso
  ncdf.variables['scama'][:,:,:,:,:,:] = scama
  # add extra variables
  print('calculate extra variables')

  ncdf.variables['sca'][:,:,:,:] = sca
  volconv = 4./3. * np.array(sizes)

  # grasp kernels need to be divided by a factor of log(x[n+1]/x[n]), i.e. the log of bin size ratios
  # for the size bins used here this ends up being a multiplication by 3.68176736
  graspfactor = 3.68176736
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

  # scama is normalized so we un-normalize it by multiplying by scattering cross-section?
  pBck = f11back / sca # get p11 from f11 by dividing by sca
  qBck = pBck * ncdf.variables['qsca'][:,:,:,:] # get backscattering efficiency by multiplying p11 by qsca

  ncdf.variables['qb'][:,:,:,:] = qBck

  angs = np.radians(angs)
  normfactor = np.sum(f11 * np.sin(angs), axis=-1)
  f11n = np.array([f11[:,:,:,:,ai] / normfactor for ai in range(len(angs))])
  g_pre = np.cos(angs) * np.sin(angs) * f11n.T
  g = np.sum(g_pre, axis=-1).T
  ncdf.variables['g'][:,:,:,:] = g

def getSizes(pfx):
  grid = 'grid1.txt'
  fn = os.path.join(pfx, grid)
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


# TODO supply path and mode from command line
fun('grasp')
