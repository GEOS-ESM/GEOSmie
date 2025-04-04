import netCDF4
import numpy as np
import bisect
from scipy.interpolate import interp1d
from scipy.stats import lognorm
import os
import numba
import hydrophobic
import sys
from pymiecoated.mie_coated import MultipleMie
import particleparams as pp

"""
These govern how integration and dimensions etc are done, so define them here
"""
scatkeys = ['s11', 's12', 's22', 's33', 's34', 's44']
#scalarkeys = ['qext', 'qsca', 'qabs', 'g', 'qb', 'csca', 'cext']
scalarkeys = ['qext', 'qsca', 'qabs', 'qb', 'g', 'csca', 'cext']
extrakeys = ['ssa', 'bsca', 'bext', 'bbck', 'refreal', 'refimag', 'lidar_ratio']
elekeys = ['pback']
nlscalarkeys = ['mass', 'volume', 'area', 'rEff', 'rMass', 'rhop', 'growth_factor', 'rLow', 'rUp']

allkeys = scatkeys + scalarkeys + extrakeys + elekeys + nlscalarkeys

globalSpheroid = None # use global variable to access spheroid data

# copy-paste from MAML_OpticsMod
def bilinear_interpolation(x1, x2, y1, y2, f1, f2, f3, f4, x, y):
  t = (x - x1) / (x2 - x1)
  u = (y - y1) / (y2 - y1)

  ret = (1.0 - t) * (1.0 - u) * f1 + \
        (      t) * (1.0 - u) * f2 + \
        (1.0 - t) * (      u) * f4 + \
        (      t) * (      u) * f3
  return ret

def get_interpolated(data, mr, mi, allmr, allmi, debug=False):
  mri = bisect.bisect_left(allmr, mr) - 1

  # clamp mri, mii 
  mri = min(len(allmr) - 2, max(0, mri))
  mii = bisect.bisect_left(allmi, mi) - 1
  mii = min(len(allmi) - 2, max(0, mii))

  if len(data.shape) == 3:
    val1 = data[mri  ,mii  ,:]
    val2 = data[mri+1,mii  ,:]
    val3 = data[mri+1,mii+1,:]
    val4 = data[mri  ,mii+1,:]
  elif len(data.shape) == 4:
    val1 = data[mri  ,mii  ,:,:]
    val2 = data[mri+1,mii  ,:,:]
    val3 = data[mri+1,mii+1,:,:]
    val4 = data[mri  ,mii+1,:,:]

  # prevent out of bounds behavior
  mr = np.clip(mr, allmr[0], allmr[-1])
  mi = np.clip(mi, allmi[0], allmi[-1])

  val     = bilinear_interpolation(allmr[mri],\
                                   allmr[mri+1],\
                                   allmi[mii],\
                                   allmi[mii+1],\
                                   val1,\
                                   val2,\
                                   val3,\
                                   val4,\
                                   mr,\
                                   mi)

  return val

def integrateShapes(data, fracs0):
  # return shape/ratio integrated data matrix
  # format is a dictionary with each value a multidimensional array
  # basically, same as the netcdf but just not a netcdf
  ret = {}
  fracs = np.array(fracs0) / np.sum(fracs0) # normalize
  for key in data.variables.keys():
    thisdata = data.variables[key][:]
    if len(thisdata.shape) == 4:
      ret[key] = np.tensordot(fracs, thisdata, axes = (0, 0))
    elif len(thisdata.shape) == 6:
      ret[key] = np.tensordot(fracs, thisdata, axes = (0, 0))
    elif len(thisdata.shape) == 1: # one-dimensional stuff, no shape averaging
      ret[key] = thisdata
    else:
      print("unhandled number of dimensions %d in key %s"%(len(thisdata.shape), key))
    
  return ret

def getDR(arr):
  diffs1 = []
  diffs1.append(arr[1] - arr[0])
  for i in range(1,len(arr)-1):
    dd1 = arr[i] - arr[i-1]
    dd2 = arr[i+1] - arr[i]
    diffs1.append((dd1+dd2) / 2.)
  diffs1.append(arr[-1] - arr[-2])
  return np.array(diffs1)

def getXArrCarma(minx, maxx, nbinperdecade):
  """
  Carma-style bins

  Best I can understand is that carma does exponential masses,
  e.g. mass = mass * x, and then gets the dr from the mass/volume

  this makes the whole decade division completely artificial and unnecessary,
  which is great as it lets us skip many of the largest size parameters and just
  use directly the min/max sizes

  so basically we would calculate decades by max/min and from that get total num of bins

  then, we'd calculate 
  rat = rMaxUse/rMinUse
  rmRat = (rat^3)^(1./double(nBinMin))
  rMin = rMinUse*((1.+rmRat)/2.)^(1.d/3)

  this means we divide the whole range into nbinmin steps (which itself is decades * binperdec)
  and apply a constant multiplier at each step, the multiplier being rmrat
  (rmass[ibin]   = rmassmin*rmrat^double(ibin))
  r[ibin]       = (rmass[ibin]/rhop/cpi)^(1.d/3.)

  main difference is that I want to use volumes and not masses so we can avoid having rho here
  thus, we get the volume of each bin (or rather the volume of the particle in that bin)
  and from that we can easily get the r using vrfact
  """

  rat = maxx / minx
  numdec = np.log10(maxx / minx) 
  nbin = numdec * nbinperdecade
  rmRat = (rat ** 3) ** (1./nbin)
  rMin = minx * ((1. + rmRat) / 2.) ** (1./3.)

  cpi = 4. / 3. * np.pi

  rvolmin = cpi * rMin ** 3.

  # convert from volume to radius bin width (dr)
  vrfact = ( (3./2./np.pi / (rmRat+1)) ** (1./3.)) * (rmRat ** (1./3.) - 1.)

  ret = []
  drarr = []
  for i in range(int(nbin)):
    rvol   = rvolmin * rmRat ** i
    r       = (rvol / cpi) ** (1./3.)
    dr      = vrfact * (rvol) ** (1./3.)
    ret.append(r)
    drarr.append(dr)

  return np.array(ret), np.array(drarr)

def getXArr(minx, maxx, nbinperdecade):
  decades = maxx / minx
  ndec = max(1,int(np.log10(decades)))
  xxrang = []
  for dd in range(ndec):
    binmin = minx * 10. ** dd
    lspace = np.linspace(binmin, 10*binmin, nbinperdecade)
    xxrang = np.concatenate([xxrang, lspace[:-1]])

  return xxrang

# return multimodal lognorm PSD
def getMMPSD(xxArr, rmodeArr, rMaxArr, rMinArr, sigmaArr, lambd, fracarr):
  totfrac = np.sum(fracarr)
  totdist = np.zeros_like(xxArr)
  for ii, rMode in enumerate(rmodeArr):
    thisdist = pp.getLogNormPSD(rMode, sigmaArr[ii], xxArr, lambd, rMaxArr[ii], rMinArr[ii])
    totdist += fracarr[ii] * thisdist
  totdist /= totfrac
  return totdist

def find_closest_ind(myList, myNumber, typ='none',ide=''):
    """
    Assumes myList is sorted. Returns closest value to myNumber.

    If two numbers are equally close, return the smallest number.
    """
    pos = bisect.bisect_left(myList, myNumber) 
    if pos == 0:
        print("WARNING! minimum index in %s"%ide)
        #sys.exit()
        return 0
    if pos == len(myList):
        print("WARNING! maximum index in %s"%ide)
        #sys.exit()
        return len(myList) - 1
    before = myList[pos - 1]
    after = myList[pos]
    # force ceil/floor behavior
    if typ == 'ceil':
      return pos
    elif typ == 'floor':
      return pos - 1

    if after - myNumber < myNumber - before:
      return pos
    else:
      return pos - 1


def createNCDF(ncdfID, oppfx, rarr, rharr, lambarr, ang, oppclassic):

  if oppclassic:
    ncdf = netCDF4.Dataset(os.path.join(oppfx, 'optics_%s.legacy.nc4'%ncdfID), 'w')
    radiusNm = 'radius'
    lambdaNm = 'lambda'
    npolNm   = 'nPol'
  else:
    ncdf = netCDF4.Dataset(os.path.join(oppfx, 'optics_%s.nc4'%ncdfID), 'w')
    radiusNm = 'bin'
    lambdaNm = 'wavelength'
    npolNm   = 'p'

  idRh     = ncdf.createDimension('rh', len(rharr))
  idLambda = ncdf.createDimension(lambdaNm, len(lambarr))
  idBin    = ncdf.createDimension(radiusNm, len(rarr))
  idNpol   = ncdf.createDimension(npolNm, 6)
  idAng    = ncdf.createDimension('ang', len(ang))

  vardict = {} # holds all the variable metadata
  vardict[radiusNm] = {'units': 'dimensionless', \
  'long_name': 'radius bin index (1-indexed)'
  }

  vardict[npolNm] = {'units': 'dimensionless', \
  'long_name': 'Scattering matrix element index, ordered as P11, P12, P33, P34, P22, P44' 
  }

  vardict['rh'] = {'units': 'fraction', \
  'long_name': 'relative humidity' 
  }
   
  vardict[lambdaNm] = {'units': 'm', \
  'long_name': 'wavelength' 
  }

  vardict['ang'] = {'units': 'degrees', \
  'long_name': 'scattering angle' 
  }

  vardict['s11'] = {'units': 'dimensionless', \
  'long_name': 'S11 element of the scattering matrix'
  }
  vardict['s12'] = {'units': 'dimensionless', \
  'long_name': 'S12 element of the scattering matrix'
  }
  vardict['s22'] = {'units': 'dimensionless', \
  'long_name': 'S22 element of the scattering matrix'
  }
  vardict['s33'] = {'units': 'dimensionless', \
  'long_name': 'S33 element of the scattering matrix'
  }
  vardict['s34'] = {'units': 'dimensionless', \
  'long_name': 'S34 element of the scattering matrix'
  }
  vardict['s44'] = {'units': 'dimensionless', \
  'long_name': 'S44 element of the scattering matrix'
  }
  vardict['qsca'] = {'units': 'dimensionless', \
  'long_name':'scattering efficiency'
  }
  vardict['qabs'] = {'units': 'dimensionless', \
  'long_name':'absorption efficiency'
  }
  vardict['qext'] = {'units': 'dimensionless', \
  'long_name':'extinction efficiency'
  }
  vardict['g'] = {'units': 'dimensionless', \
  'long_name': 'asymmetry factor'
  }
  vardict['ssa'] = {'units': 'dimensionless', \
  'long_name': 'single-scattering albedo'
  }
  vardict['bsca'] = {'units': 'm2 (kg dry mass)-1', \
  'long_name': 'mass scattering efficiency'
  }
  vardict['bext'] = {'units': 'm2 (kg dry mass)-1', \
  'long_name': 'mass extinction efficiency'
  }
  vardict['csca'] = {'units': 'm2', \
  'long_name': 'mass scattering cross-section'
  }
  vardict['cext'] = {'units': 'm2', \
  'long_name': 'mass extinction cross-section'
  }
  vardict['bbck'] = {'units': 'm2 (kg dry mass)-1', \
  'long_name': 'mass backscatter efficiency'
  }
  vardict['lidar_ratio'] = {'units': '', \
  'long_name': 'lidar ratio'
  }
  vardict['mass'] = {'units': 'kg', \
  'long_name': 'particle mass'
  }
  vardict['volume'] = {'units': 'm3 kg-1', \
  'long_name': 'particle volume per kg dry mass'
  }
  vardict['area'] = {'units': 'm2 kg-1', \
  'long_name': 'particle cross sectional area per kg dry mass'
  }
  vardict['rEff'] = {'units': 'm', \
  'long_name': 'effective radius of bin'
  }
  vardict['rMass'] = {'units': 'kg', \
  'long_name': 'effective mass of wet particle'
  }
  vardict['rUp'] = {'units': 'm', \
  'long_name': 'upper edge of radius bin'
  }
  vardict['rLow'] = {'units': 'm', \
  'long_name': 'lower edge of radius bin'
  }

  vardict['pback'] = {'units': 'dimensionless', \
  'long_name': 'phase function in backscatter direction, ordered as P11, P12, P33, P34, P22, P44' 
  }
  vardict['qb'] = {'units': 'dimensionless', \
  'long_name': 'backscattering efficiency' 
  }
  vardict['rhop'] = {'units': 'kg m-3', \
  'long_name': 'wet particle density'
  }
  vardict['growth_factor'] = {'units': 'fraction', \
  'long_name': 'growth factor = ratio of wet to dry particle radius'
  }
  vardict['refreal'] = {'units': 'dimensionless', \
  'long_name': 'real refractive index of wet particle'
  }
  vardict['refimag'] = {'units': 'dimensionless', \
  'long_name': 'imaginary refractive index of wet particle'
  }

  idRH     = ncdf.createVariable('rh', 'f8', ('rh'), compression='zlib')
  idLambda = ncdf.createVariable(lambdaNm, 'f8', (lambdaNm), compression='zlib')
  idR      = ncdf.createVariable(radiusNm, 'i8', (radiusNm), compression='zlib')
  idNpol   = ncdf.createVariable(npolNm, 'i8', (npolNm), compression='zlib')
  idAng    = ncdf.createVariable('ang', 'f8', ('ang'), compression='zlib')

  # use list for easier looping
  dimvars = [radiusNm, 'rh', lambdaNm, 'ang', npolNm]
  for var in dimvars:
    ncdf.variables[var].long_name = vardict[var]['long_name']
    ncdf.variables[var].units = vardict[var]['units']

  # idR is a 1-indexed bin index
  idR[:] = range(1,len(rarr)+1)
  idRH[:] = rharr
  idLambda[:] = lambarr
  idAng[:] = ang
  idNpol[:] = [11,12,33,34,22,44]

  if oppclassic:
    ncdfDimensionTypes = {\
                          "scalarScattering": (radiusNm, "rh", lambdaNm),\
                          "noLambdaScalarScattering": (radiusNm, "rh"),\
                          "scatteringFun": (radiusNm, "rh", lambdaNm, "ang"),\
                          "elementScattering": (npolNm, radiusNm, "rh", lambdaNm),\
                          "scalarVals": (radiusNm),\
                          }
  else:
    ncdfDimensionTypes = {\
                          "scalarScattering": (radiusNm, lambdaNm, "rh"),\
                          "noLambdaScalarScattering": (radiusNm, "rh"),\
                          "scatteringFun": (radiusNm, lambdaNm, "rh", "ang"),\
                          "elementScattering": (radiusNm, lambdaNm, "rh", "p"),\
                          "scalarVals": (radiusNm),\
                          }

  # for each NetCDF variable, assign one dimension type based on the dimension types dictionary
  ncdfVariableDimensionTypes = {\
  "s11": "scatteringFun",\
  "s12": "scatteringFun",\
  "s22": "scatteringFun",\
  "s33": "scatteringFun",\
  "s34": "scatteringFun",\
  "s44": "scatteringFun",\
  "qsca": "scalarScattering",\
  "qabs": "scalarScattering",\
  "qext": "scalarScattering",\
  "g": "scalarScattering",\
  "ssa": "scalarScattering",\
  "qb": "scalarScattering",\
  "bsca": "scalarScattering",\
  "bext": "scalarScattering",\
  "csca": "scalarScattering",\
  "cext": "scalarScattering",\
  "bbck": "scalarScattering",\
  "lidar_ratio": "scalarScattering",\
  "mass": "noLambdaScalarScattering",\
  "volume": "noLambdaScalarScattering",\
  "area": "noLambdaScalarScattering",\
  "rEff": "noLambdaScalarScattering",\
  "rMass": "noLambdaScalarScattering",\
  "rUp": "noLambdaScalarScattering",\
  "rLow": "noLambdaScalarScattering",\
  "pback": "elementScattering",\
  "rhop": "noLambdaScalarScattering",\
  "growth_factor": "noLambdaScalarScattering",\
  "refreal": "scalarScattering",\
  "refimag": "scalarScattering",\
  }


  ncdfVars = {}
  ncdfVarDims = {}

  # assign dimensions to each NetCDF variable from the dimension type dictionary
  for ncdfKey, dimensionKey in list(ncdfVariableDimensionTypes.items()):
    ncdfVarDims[ncdfKey] = ncdfDimensionTypes[dimensionKey]

  # create NetCDF variables 
  # TODO define the variable types somewhere instead of using f8 for all
  for ncdfKey in list(ncdfVariableDimensionTypes.keys()):
    ncdfVars[ncdfKey] = ncdf.createVariable(ncdfKey, 'f8', ncdfVarDims[ncdfKey], compression='zlib')
    ncdf.variables[ncdfKey].long_name = vardict[ncdfKey]['long_name']
    ncdf.variables[ncdfKey].units = vardict[ncdfKey]['units']

  return ncdf


def initializeXarr(params, radind, minlam, maxlam):
  ###
  ### Set up xxarr and drarr 
  ###
  pparam = params['psd']['params']
  psdtype = params['psd']['type']
  if psdtype == 'lognorm':
    # we should also consider max-humidity case for maxx0

    # use minx and maxx instead of weird r0 system
    # PRC: Not sure why OPK made limmaxx0 greater by factor of 3
    limminx0 = pparam['rmin0'][radind][0] * 2 * np.pi / maxlam
    limmaxx0 = pparam['rmax0'][radind][-1] * 2 * np.pi / minlam * 3.0

    if 'numperdec' in pparam:
      numperdec = pparam['numperdec'][radind]
    else:
      print("numperdec missing from dict")
      sys.exit()

    xxarr, drarr = getXArrCarma(limminx0, limmaxx0, numperdec)
    drarr2 = getDR(xxarr)
    drrat = drarr / drarr2

  elif psdtype == 'ss':
    minx = pparam['rMinMaj'][radind] * 2 * np.pi / maxlam 
    maxx = pparam['rMaxMaj'][radind] * 2 * np.pi / minlam * 10 # account for relative humidity growth?
    numperdec = pparam['numperdec'][radind]
    xxarr, drarr = getXArrCarma(minx, maxx, numperdec) 
    drarr2 = getDR(xxarr)
    drrat = drarr / drarr2
  elif psdtype == 'du':
    minx = pparam['rMinMaj'][radind][0] * 2 * np.pi / maxlam
    maxx = pparam['rMaxMaj'][radind][-1] * 2 * np.pi / minlam
    # Try a logarithmic spacing
    xxarr = (np.linspace(np.log10(minx), np.log10(maxx), 1000))**10.
    drarr = None # not used for dust since it's overwritten by GRASP values
    
  return xxarr, drarr
    
def copyDryValues(opncdf, allkeys, scatkeys, elekeys, extrakeys, nlscalarkeys, radind, rhi, li, oppclassic):
  if oppclassic:
    for key in allkeys: 
      if key in scatkeys:
        opncdf.variables[key][radind, rhi, li, :] = opncdf.variables[key][radind, 0, li, :]
      elif key in elekeys:
        opncdf.variables[key][:, radind, rhi, li] = opncdf.variables[key][:, radind, 0, li]
      elif key in scalarkeys + extrakeys:
        opncdf.variables[key][radind, rhi, li] = opncdf.variables[key][radind, 0, li]
      elif key in nlscalarkeys:
        opncdf.variables[key][radind, rhi] = opncdf.variables[key][radind, 0]
      else:
        print("key category missing: %s"%key)
        sys.exit()
  else:
    for key in allkeys: 
      if key in scatkeys:
        opncdf.variables[key][radind, li, rhi, :] = opncdf.variables[key][radind, li, 0, :]
      elif key in elekeys:
        opncdf.variables[key][radind, li, rhi, :] = opncdf.variables[key][radind, li, 0, :]
      elif key in scalarkeys + extrakeys:
        opncdf.variables[key][radind, li, rhi] = opncdf.variables[key][radind, li, 0]
      elif key in nlscalarkeys:
        opncdf.variables[key][radind, rhi] = opncdf.variables[key][radind, 0]
      else:
        print("key category missing: %s"%key)
        sys.exit()

def getHumidRefractiveIndex(params, radind, rhi, rh, nref0, nrefwater):
  psdtype = params['psd']['type']
  pparam = params['psd']['params']
  rparams = params['rhDep']
  if params['rhDep']['type'] == 'simple' or (params['rhDep']['type'] == 'trivial' and rhi == 0):
    rhparams = params['rhDep']['params']
    simplegf = rhparams['gf'][rhi]

    gf = simplegf
    rrat = (1./simplegf)

    # formula from Pete's runmie.pro
    nrefUse = [nrefwater + (nref0[i] - nrefwater) * (rrat) ** 3. for i in range(len(nref0))]
    mr = [nrefUse[i].real for i in range(len(nrefUse))] 
    mi = [nrefUse[i].imag for i in range(len(nrefUse))]
  elif params['rhDep']['type'] == 'ss':
    rMinMaj = pparam['rMinMaj'][radind] 
    onerh = rh[rhi]
    rMinUse = pp.humidityGrowth(rparams, rMinMaj, onerh, rh)
    rrat = (rMinMaj/rMinUse) # ratio of dry binmin and wet binmin
    gf   = 1./rrat # inverse of the ratio, i.e. linear growth factor
    nrefUse = [nrefwater + (nref0[i] - nrefwater) * (rrat) ** 3. for i in range(len(nref0))]
    mr = [nrefUse[i].real for i in range(len(nrefUse))] 
    mi = [nrefUse[i].imag for i in range(len(nrefUse))]
    
  else:
    print("PROBLEM! NO RHDEP DEFINED!")
    sys.exit()

  return mr, mi, gf, rrat


def calculatePSD(params, radind, onerh, rh, xxarr, drarr, rrat, lam):
  pparam = params['psd']['params']
  psdtype = params['psd']['type']
  rparams = params['rhDep']

  xconv = 2 * np.pi / lam
  rarr = xxarr / xconv

  # i'm not super happy about doing this with this sort of if-else
  # consider refactoring to a cleaner solution
  if (psdtype == 'lognorm'):
    rmodes = [pp.humidityGrowth(rparams, oner0, onerh, rh) for oner0 in pparam['r0'][radind]]
    rmaxs0 =  [onermax0 for onermax0 in pparam['rmax0'][radind]]
    rmaxs =  [pp.humidityGrowth(rparams, onermax0, onerh, rh) for onermax0 in rmaxs0]
    # rmin does not grow with humidity intentionally
    rmins =  [onermin0 for onermin0 in pparam['rmin0'][radind]]

    rMinUse = rmins[0]
    rMaxUse = rmaxs[0]
    rLow = rmins[0]
    rUp = rmaxs0[0]

    # no humidity growth for sigma
    sigmas = pparam['sigma'][radind] 

    fracs = pparam['fracs']
    psd = []
    ref = []
    for psdi in range(len(rmodes)):
#     This call returns PSD as dN/dr
      thispsd = pp.getLogNormPSD(rmodes[psdi], sigmas[psdi], xxarr, lam, rmaxs[psdi], rmins[psdi])
      psdbins = drarr

      dr = psdbins
      dndr = thispsd / np.sum(thispsd)

#     Following converts PSD to dN      
      thispsd *= psdbins 
      thispsd /= np.sum(thispsd)

      psd.append(thispsd)

      rrarr = xxarr * lam / 2. / np.pi
      thisref = np.sum(rrarr**4. * thispsd) / np.sum(rrarr**3. * thispsd)
      ref.append(thisref)

  elif psdtype == 'ss':
    rMinMaj = pparam['rMinMaj'][radind] 
    rMaxMaj = pparam['rMaxMaj'][radind]
    rLow = rMinMaj
    rUp = rMaxMaj
    rMinUse = pp.humidityGrowth(rparams, rMinMaj, onerh, rh)
    rMaxUse = pp.humidityGrowth(rparams, rMaxMaj, onerh, rh)

    nBinMin = len(rarr)
    rrat = (rMinMaj/rMinUse) # ratio of dry binmin and wet binmin
    rmRat = (rrat ** 3) ** (1. / nBinMin)

    small = np.where(rarr < rMinUse)[0]
    large = np.where(rarr > rMaxUse)[0]


    psdbins = getDR(rarr)
    dr = psdbins

    r80rat = 1.65 * rrat      # ratio of the r80 radius to the wet radius
    r80  = rarr * r80rat * 1e6  # radius in r80 space in um
    dr80 = dr * r80rat * 1e6

    aFac = 4.7 * (1. + 30. * r80) ** (-0.017 * r80 ** (-1.44))
    bFac = (0.433 - np.log10(r80)) / 0.433
    dndr80 = 1.373 * r80 ** (-aFac) * (1. + 0.057 * r80 ** 3.45) * 10. ** (1.607 * np.exp(-bFac ** 2.))
    dndr = dndr80 * r80rat 

    # enforce bin limits
    dndr[small] = 0.
    dndr[large] = 0.

    # normalize dndr 
    dndr = dndr / np.sum(dndr)

    psd = np.copy(dndr) * dr
    psd /= np.sum(psd) 
    psd = [psd] # force into array to match multimodal PSD system

    ref = [np.sum((xxarr * lam / 2. / np.pi)**4 * psd) / np.sum((xxarr * lam / 2. / np.pi)**3 * psd)]

  elif psdtype == 'du':
    psd = []
    ref = []
    numminor = len(pparam['rMinMaj'][radind])
    for rMinorInd in range(numminor):
      rMinMaj = pparam['rMinMaj'][radind][rMinorInd] 
      rMaxMaj = pparam['rMaxMaj'][radind][rMinorInd]
      rMinUse = rMinMaj
      rMaxUse = rMaxMaj
      rLow = rMinMaj
      rUp = rMaxMaj
      small = np.where(rarr < rMinMaj)[0]
      large = np.where(rarr > rMaxMaj)[0]

      dndr = rarr ** (-4.)
      psdbins = getDR(xxarr) 
  
      dndr[small] = 0.
      dndr[large] = 0.

      dr = psdbins

      if(np.max(dndr) == 0.):
        print(rarr)
        print(xxarr)
        print(lam)
        print('prc:',np.min(xxarr),np.max(xxarr),np.min(rarr),np.max(rarr),lam,rMinMaj,rMaxMaj)
        sys.exit()


      thispsd = np.copy(dndr) * dr
      thispsd /= np.sum(thispsd)
      psd.append(thispsd)

      rrarr = xxarr * lam / 2. / np.pi
      thisref = np.sum(rrarr**4. * thispsd) / np.sum(rrarr**3. * thispsd)
      ref.append(thisref)

  return psd, ref, rLow, rUp

# This is the main integration function called from runoptics.py
# Inputs are:
#  partID0:    the particle type JSON file header or filename, e.g., bc, oc, ...
#  datatype:   the type of the file to be parsed, must be JSON presently
#  oppfx:      the path for the output file
#  oppclassic: generate legacy format lookup table
def fun(partID0, datatype, oppfx, oppclassic):

  # clean up partID, reduce just to particle type (e.g., bc, oc, ...)
  partID = partID0.split('/')[-1].replace(".json", "")

  print("\n ####################\n Starting case %s\n ####################\n"%partID)

  ncdfID = '%s'%partID
  if '-orig' in partID:
    partID2 = partID0.replace('-orig', '')
  else:
    partID2 = partID0
  # Get the particle properties from the JSON file
  params = pp.getParticleParams(partID2, datatype)

  # Particle shape determines calculation path, either Mie calculations
  # or using GRASP-like kernel files. Code does not presently include
  # any alternative internal calculations to Mie.
  mode = 'mie'
  if 'mode' in params:
    mode = params['mode']
#  if mode == 'spheroid' or mode == 'spheroid_sphere':
  if mode == 'kernel':
    useGrasp = True
  elif mode == 'mie':
    useGrasp = False
    
  # Refractive indices for particles and water
  mList = params['mList']
  waterMList = pp.getWaterM()

  # List of wavelengths in the particle refractive index table
  allLambda = mList[0][0]

  # interpolation functions for refractive indices mr, mi
  partMr = [interp1d(mList[i][0], mList[i][1]) for i in range(len(mList))]
  partMi = [interp1d(mList[i][0], mList[i][2]) for i in range(len(mList))]
  waterMr = interp1d(waterMList[0], waterMList[1])
  waterMi = interp1d(waterMList[0], waterMList[2])


  radindarr = None
  radiusarr = None
  psdtype = params['psd']['type']
  if psdtype == 'lognorm':
    radindarr = list(range(len(params['psd']['params']['r0'])))
    radiusarr = list(params['psd']['params']['r0'])
  elif psdtype == 'ss':
    radindarr = list(range(len(params['psd']['params']['rMinMaj'])))
    radiusarr = list(params['psd']['params']['rMinMaj'])
  elif psdtype == 'du':
    radindarr = list(range(len(params['psd']['params']['rMinMaj'])))
    radiusarr = [xx[0] for xx in params['psd']['params']['rMinMaj']]


  # Wavelengths to compute on, presently defined by set of wavelengths
  # defined in the particle refractive index files
  lambarr = allLambda

  rh = params['rh'] 

  """
  Define parameters for netcdf creation
  """

  # Angles phase functions will be written at
  if mode =='mie':
    # Define output scattering angles
    ang1 = np.linspace(0., 1., 100, endpoint=False)
    ang2 = np.linspace(1., 10., 100, endpoint=False)
    ang3 = np.linspace(10., 180., 171, endpoint=True)
    ang = np.concatenate([ang1,ang2,ang3])
  elif useGrasp:
    ang = np.linspace(0., 180., 181)

  costarr = np.cos(np.radians(ang))
  minlam = lambarr[0]
  maxlam = lambarr[-1]

  opncdf = createNCDF(ncdfID, oppfx, radiusarr, rh, lambarr, ang[:], oppclassic)

  if useGrasp:
    # if we are using a spheroid kernel system then override xxarr with
    # what is actually available from the spheroids

    kparams = params['kernel_params']

    if 'path' not in kparams:
      print('kernel path parameter (\'path\') not defined')
      sys.exit()
    if 'shape_dist' not in kparams:
      print('kernel shape distribution parameter (\'shape_dist\') not defined')
      sys.exit()

    spdata = readSpheroid(params['kernel_params']['path'])
    distpath = kparams['shape_dist']
    spfracs = np.loadtxt(distpath, usecols=[0], unpack=True, ndmin=1)
    print('Integrating kernels...')
    globalSpheroid = integrateShapes(spdata, spfracs)
    print('Done')

    xxarr = spdata.variables['x'][:]
    # Original call to getDR below, but for our kernel files to date
    # there is a simple geometric progression in size space; i.e., 
    # x1 = rat*x0, x2 = rat*x1, ...
    # and so more accurately r/dr is constant. In CARMA land I reproduce
    # some of the functionality here to get this accurate.
    # drarr = getDR(xxarr) 
    rmrat = (xxarr[1]/xxarr[0])**3
    vrfact = ( (3./2./np.pi / (rmrat+1))**(1./3.))*(rmrat**(1./3.) - 1.)
    drarr  = vrfact*(4./3.*np.pi*xxarr**3.)**(1./3.)

  # Loop over the particle size bins/modes
  for radind in radindarr:
    print("=== === === USING RADIND %d"%radind)

    if not useGrasp:
      xxarr, drarr = initializeXarr(params, radind, minlam, maxlam)
    else:
      # Get a nominal size array like you are not using GRASP for later
      xxarr_, drarr_ = initializeXarr(params, radind, minlam, maxlam)

    if mode == 'mie':
      multipleMie = MultipleMie(xxarr, None, costarr)
      multipleMie.preCalculate()

    allvals = {}
    for key in allkeys:
      allvals[key] = np.zeros(opncdf.variables[key][:].shape)

    """
    Start wavelength loop
    TODO!
    parallelization over lambda, i.e. have a single worker evaluate each lambda since they are independent of each other
    therefore, we should move this huge block of code under the loop into a separate function that takes lambda as a 
    parameter along with everything else it needs
    """
    for li, lam in enumerate(lambarr):
      print("+++++ LAMBDA %.2e +++++"%lam)
      mr0 = [partMr[i](lam) for i in range(len(partMr))]
      mi0 = [-partMi[i](lam) for i in range(len(partMi))]
      #mi0 = -partMi(lam) # defined as negative 
      nref0 = [complex(mr0[i], mi0[i]) for i in range(len(mr0))]

      xconv = 2 * np.pi / lam
      rarr = xxarr / xconv

      if params['rhDep']['type'] == 'trivial':
        nrefwater = complex(1,0) # placeholder, never used
      else:
        watermr0 = waterMr(lam)
        watermi0 = waterMi(lam)
        nrefwater = complex(watermr0, watermi0)

      if 'maxrh' in params:
        maxrh = params['maxrh']
        capind = np.where(np.array(rh) > maxrh)[0]
        rh = np.array(rh)
        rh[capind] = maxrh

      """
      Get Dry Particle Properties
      """
      mr0, mi0, gf, rrat0 = getHumidRefractiveIndex(params, radind, 0, rh, nref0, nrefwater)
      psd0, reff_mass0, rLow0, rUp0 = calculatePSD(params, radind, 0., rh, xxarr, drarr, rrat0, lam)


      """
      Start RH loop
      """
      for rhi, onerh in enumerate(rh):
        pparam = params['psd']['params']
        rparams = params['rhDep']

        if params['rhDep']['type'] == 'trivial' and rhi > 0: # same values for all rh, save in computation
          copyDryValues(opncdf, allkeys, scatkeys, elekeys, extrakeys, nlscalarkeys, radind, rhi, li, oppclassic)
          continue

        mr, mi, gf, rrat = getHumidRefractiveIndex(params, radind, rhi, rh, nref0, nrefwater)

        psd, ref, rLow, rUp = calculatePSD(params, radind, onerh, rh, xxarr, drarr, rrat, lam)
        if useGrasp:
          psd_, ref_, rLow_, rUp_ = calculatePSD(params, radind, onerh, rh, xxarr_, drarr_, rrat, lam)

        """
        ***********
        
        Start the calculations

        ***********
        """

        rhop00 = params['rhop0'] # read from a file
        if isinstance(rhop00, list): 
          # rhop0 is defined separately for each size bin, read the right one
          rhop0 = rhop00[radind]
        else:
          rhop0 = rhop00

        rhow  = 1000. # density of water, constant
        rhop = rrat ** 3. * rhop0 + (1. - rrat ** 3.) * rhow

        if mode == 'mie':
          allret = []
          for refi in range(len(mr)):
            rawret = rawMie(multipleMie, scatkeys, scalarkeys, lam, mr[refi], mi[refi], None, costarr)
            allret.append(rawret)
          
          if len(allret) == 1:
            # make compatible with multibin psd
            allret = [allret[0] for i in range(len(psd))] # multibin

          retkeys = list(allret[0].keys()) # all are assumed to have the identical keys so we just get them from 0th index

          # separate integration step
          ret = integratePSD(multipleMie.xArr, allret, psd, pparam['fracs'][radind], lam, reff_mass0, rhop0, rhop)

        elif useGrasp:
          ret0 = globalSpheroid

          allmr = ret0['mr'][:]
          allmi = -ret0['mi'][:] # change sign

          keys = ['ext', 'abs', 'sca', 'qext', 'qabs', 'qsca', 'qb', 'g', 'cext', 'csca', 'cabs']
#          keys = ['ext', 'abs', 'sca', 'qext', 'qsca', 'qb', 'g', 'cext', 'csca']
          scatelekeys = ['scama']
          scatelems = ['s11', 's22', 's33', 's44', 's12', 's34'] # order Mischenko's code expects
          ret1 = {}
          for kk in keys:
            valmat = ret0[kk][:,:,:]
            val = get_interpolated(valmat, mr, mi, allmr, allmi)
            ret1[kk] = val

          for kk in scatelekeys:
            for scati in range(len(scatelems)):
              if scati == 0:
                dodebug = True
              else:
                dodebug = False
              valmat = ret0[kk][:,:,:,scati,:]
              val = get_interpolated(valmat, mr, mi, allmr, allmi, debug=dodebug)
              ret1[scatelems[scati]] = val

          ret1['qsca'] = ret1['qsca']

          ret1['csca'] = np.array(ret1['qsca']) * np.pi * rarr ** 2
          ret1['cext'] = np.array(ret1['qext']) * np.pi * rarr ** 2

          ret1 = [ret1 for i in range(len(psd))] # multibin
          if len(pparam['fracs']) == 1:
            # if only one set of fracs is given we always use it
            fracs = pparam['fracs'][0]
          else:
            fracs = pparam['fracs'][radind]
          retkeys = list(ret1[0].keys()) # all are assumed to have the identical keys so we just get them from 0th index
          ret = integratePSD(xxarr, ret1, psd, fracs, lam, reff_mass0, rhop0, rhop)

          # This is a hack and only applied if (a) psd type is lognormal or 
          # (b) psd type is 'du' and the number fraction = 1 (i.e., not the first bin)
          # Post integration rescaling of the mass efficiencies because of limited resolution of kernel tables
          # introducing error in effective radius calculation
          # Note: this could be done more generally for number, volume, etc that should not vary with wavelength
          # For now keep it simple and fix the extinction efficiencies
          if(psdtype == 'lognorm') or (psdtype == 'du' and len(pparam['fracs'][radind]) == 1):
            rrarr = xxarr_ * lam / (2. * np.pi)
            for fraci, frac in enumerate(fracs):
              rarr2 = rrarr ** 2. * psd_[fraci]
              rarr3 = rrarr ** 3. * psd_[fraci]
              reff = np.sum(rarr3) / np.sum(rarr2)

              ret['bsca'] = ret['bsca'] * ret['rEff']/reff
              ret['bext'] = ret['bext'] * ret['rEff']/reff
              ret['bbck'] = ret['bbck'] * ret['rEff']/reff
              ret['rEff'] = reff


        qsca = np.array(ret['qsca'])
        qext = np.array(ret['qext'])
        qb = np.array(ret['qb'])
        g  = np.array(ret['g'])

        """
        Calculate extra post-integration variables
        """

        ret['lidar_ratio'] = qext / qb * 4 * np.pi

        ret['ssa'] = qsca / qext

        ret['rLow'] = rLow
        ret['rUp'] = rUp

        pbackorder = ['s11', 's12', 's33', 's34', 's22', 's44']
        ret['pback'] = np.array([ret[key][-1] for key in pbackorder])

        ret['growth_factor'] = gf

        mass0 = ret['volume'] * rhop0

        ret['rhop'] = rhop
        ret['area'] = ret['area'] / mass0
        ret['volume'] = ret['volume'] / mass0

        # mr and mi are arrays since we might have multimodal PSD with different components as modes
        ret['refreal'] = mr[0]
        ret['refimag'] = -np.abs(mi[0]) # force negative to be consistent with Pete's tables

        pback = np.array(ret['pback'])

#       Support for legacy format lookup tables
        if oppclassic:
          for key in allkeys: # we can also save to allvals and write later
            if key in scatkeys:
              opncdf.variables[key][radind, rhi, li, :] = ret[key][:]
            elif key in elekeys:
              opncdf.variables[key][:, radind, rhi, li] = ret[key][:]
            elif key in scalarkeys + extrakeys:
              opncdf.variables[key][radind, rhi, li] = ret[key]
            elif key in nlscalarkeys:
              opncdf.variables[key][radind, rhi] = ret[key]
            else:
              print("key category missing: %s"%key)
              sys.exit()
        else:
          for key in allkeys: # we can also save to allvals and write later
            if key in scatkeys:
              opncdf.variables[key][radind, li, rhi, :] = ret[key][:]
            elif key in elekeys:
              opncdf.variables[key][radind, li, rhi, :] = ret[key][:]
            elif key in scalarkeys + extrakeys:
              opncdf.variables[key][radind, li, rhi] = ret[key]
            elif key in nlscalarkeys:
              opncdf.variables[key][radind, rhi] = ret[key]
            else:
              print("key category missing: %s"%key)
              sys.exit()
        # end rh loop
      # end lambda loop
    # end radind loop

  opncdf.close()

"""
Run Mie directly at the desired values

Calculates the Mueller matrix values from S12
"""

def calculateScatVals(vals, scatkeys, ret, numsizes, s1arr, s2arr):
  ret[0,:,:] = 0.5 * ( np.abs(s1arr) ** 2 + np.abs(s2arr) ** 2 )  # s11
  ret[1,:,:] = -0.5 * ( -np.abs(s1arr) ** 2 + np.abs(s2arr) ** 2 ) # s12
  ret[2,:,:] = 0.5 * ( np.abs(s1arr) ** 2 + np.abs(s2arr) ** 2 ) # s22 - same as s11
  ret[3,:,:] = ( s1arr * np.conj(s2arr) ).real # s33
  ret[4,:,:] = -( np.conj(s1arr) * s2arr ).imag # s34
  ret[5,:,:] = ( s1arr * np.conj(s2arr) ).real # s44 - same as s33

def readSpheroid(fn):
  print('Start reading kernels')
  if not os.path.exists(fn):
    sys.exit("Kernel path %s does not exist"%fn)
  data = netCDF4.Dataset(fn, 'r')
  print('Done')
  return data

"""
For multimodal PSD both rawret as well as psd are lists, so in the integration phase
we iterate over them
"""
def integratePSD(xxarr, rawret, psd, fracs, lam, reff0, rhop0, rhop):
  ret = {}
  rrarr = xxarr * lam / (2. * np.pi)

  ret['num'] = 0.
  ret['area'] = 0.
  ret['volume'] = 0.
  ret['mass'] = 0.
  ret['rEff'] = 0.
  ret['rMass'] = 0.
  ret['bsca'] = 0.
  ret['bext'] = 0.
  ret['bbck'] = 0.
  ret['lidar_ratio'] = 0.

  totarea = 0

  retkeys = list(rawret[0].keys()) # all are assumed to have the identical keys so we just get them from 0th index

  # initialize ret with all keys needed
  for key in retkeys:
    if key in scatkeys:
      shape = rawret[0][key].shape
      if len(shape) > 1:
        ret[key] = np.zeros(shape[1])
      else:
        ret[key] = np.zeros(shape[0])
    else:
      ret[key] = 0.

  for fraci, frac in enumerate(fracs):
    thisret = {}
    for key in ret.keys():
      thisret[key] = np.zeros_like(ret[key])

    drdndr = psd[fraci]
    reff_mass0 = reff0[fraci]

    num = np.sum(drdndr)
    rarr2 = rrarr ** 2. * drdndr
    rarr3 = rrarr ** 3. * drdndr
    rarr4 = rrarr ** 4. * drdndr

    area = np.pi * np.sum(rarr2) / num 
    volume = 4./3. * np.pi * np.sum(rarr3) / num 
    mass = volume * rhop

    thisarea = np.pi * rrarr ** 2.
    thisweight = thisarea # scale by area

    reff = np.sum(rarr3) / np.sum(rarr2)

    thisret['num'] = num
    thisret['area'] = area
    thisret['mass'] = mass
    thisret['volume'] = volume
    thisret['rEff'] = reff

    reff_mass = np.sum(rarr4) / np.sum(rarr3)
    thisret['rMass'] = 4. / 3. * np.pi * rhop * reff_mass ** 3.
    rMass0 = 4. / 3. * np.pi * rhop0 * reff_mass0 ** 3

    for key in retkeys:

#     PRC: qb returned from the mie call is the backscatter efficiency 
#     (to a factor of 2) consistent with what gets calculated from MIEV
#     but Osku's original integration is inconsistent (and wrong?) with
#     my IDL code. Since qb = p11*qsca/(4pi) I'm right here backing out
#     P11 at each subbin, which I will integrate separately
      qbwght = ['qb']
      if key in qbwght:
        thisq    = np.array(rawret[fraci]['qsca'])
        beta     = np.sum(thisq*rarr2)*np.pi
        k        = 2.*np.pi/lam
        k2       = k*k
        fac      = 4.*np.pi/(2.*k2*beta)
        fac      = 1./beta
        thisqb   = np.array(rawret[fraci][key])
#        thisp11  = thisqb/(thisq)#*(4.*np.pi)
#        thisp11  = thisqb*thisarea/np.pi
        thisp11  = thisqb*thisarea*(4.*np.pi)
        p11back  = np.dot(thisp11.T,psd[fraci])*fac#/np.sum(rarr2)
#        print(p11back,beta,np.sum(psd[fraci])*fac)
#        print(thisq)
#        print(rrarr*1e6)
#        print(thisp11)
#        exit()

#     PRC: The following seems dangerous as it resets "thisweight" and
#     seems it would be bad if the order of the keys changed
      qscawe = ['g']

      if key in qscawe:
        thisweight *= rawret[fraci]['qsca'] # scale also by qext

#     PRC: this seems to be consisent with how "thisarea" is modified above
      sumarea = np.sum(psd[fraci] * thisweight)
#     PRC: WTF, this accumulates over the loop of keys and seems totally wrong!
      totarea += frac * sumarea
      
      if key in scatkeys:
        thisvals = rawret[fraci][key]
        thisret[key] = np.dot(thisvals.T, psd[fraci]) # integrate
      elif key in qbwght:
        thisret[key] = p11back*thisret['qsca']/(4.*np.pi)
#        print(thisret[key])
#        exit()
      else:
        """
        For external mixing we first need to integrate over the individual PSD, 
        then sum the two cases 
        """

        # mixing mode has to be currently changed manually here
        #mixing = 'internal'
        mixing = 'dumix'
        #mixing = 'external'
        thisvals = np.array(rawret[fraci][key])

        if mixing == 'internal':
          thisret[key] = np.dot(thisvals, psd[fraci] * thisweight)
        elif mixing == 'external':
          thisint = np.dot(thisvals, psd[fraci] * thisweight)
          thisarea = np.sum(psd[fraci] * thisweight)
          thisret[key] = thisint 
        elif mixing == 'dumix':
          thisint = np.dot(thisvals, psd[fraci] * thisweight)
          thisret[key] = thisint / sumarea

    # calculate mass efficiencies here
    massConversion = 1. / rhop / thisret['rEff'] * thisret['rMass'] / rMass0
    thisret['bsca'] = 3./4. * thisret['qsca'] * massConversion
    thisret['bext'] = 3./4. * thisret['qext'] * massConversion
    thisret['bbck'] = 3./4. * thisret['qb'] * massConversion / (4*np.pi)

    # save all from thisret to ret
    for key in thisret.keys():
      ret[key] += frac * thisret[key]

  # post integral normalization
  # only normalize retkeys, not everything in ret
  # e.g. rEff we don't normalize since it's calculated separately
  for key in retkeys + ['bsca', 'bext', 'bbck']:
   if mixing == 'internal': 
     ret[key] /= totarea

  return ret

def rawMie(mm, scatkeys, scalarkeys, lam, mr, mi, psd, costarr):
  ret = {}
  xxarr = mm.xArr
  rrarr = xxarr[:] * lam / (2. * np.pi)

  vals0 = mm.calculateS12SizeRange(mr, mi)

  qsca = vals0['qsca']
  vals0['csca'] = np.array(vals0['qsca']) * np.pi * rrarr ** 2
  vals0['cext'] = np.array(vals0['qext']) * np.pi * rrarr ** 2
  vals0['g'] = np.array(vals0['asy']) # change name
  
  numsizes = len(vals0['s12'])
  numang = len(vals0['s12'][0])
  retvals = np.empty([6, numsizes, numang], dtype=np.float64)

  s1arr = np.zeros([numsizes, numang], dtype=np.complex128)
  s2arr = np.zeros([numsizes, numang], dtype=np.complex128)
  for sizi, size in enumerate(range(numsizes)):
    s1 = np.array([foo[0] for foo in vals0['s12'][sizi]])
    s2 = np.array([foo[1] for foo in vals0['s12'][sizi]])
    s1arr[sizi, :] = s1
    s2arr[sizi, :] = s2

  calculateScatVals(vals0['s12'], scatkeys, retvals, len(s1arr), s1arr, s2arr)

  k = 2 * np.pi / lam

  vals = {}
  for key in scalarkeys:
    vals[key] = vals0[key]
  for ki, key in enumerate(scatkeys):
    vals[key] = retvals[ki]

  keys = scatkeys + scalarkeys
  for key in keys:
    if key in scatkeys:
      thisvals = vals[key]# / bsca
    else:
      thisvals = vals[key]

    ret[key] = thisvals

  return ret
