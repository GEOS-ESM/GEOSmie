"""
 Utility functions that return particle property definitions
 defined in the JSON configuration file in a simple data
 structure. Included are simple functions to read wavelength
 dependent refractive index tables, compute humidified 
 particle sizes, and compute simple size distribution
 functions (e.g., lognormal, could consider adding gamma
 function as a possible add on.
"""

import numpy as np
import os
import json
import sys
import carma_utils

# Load and return data structure with contents of input JSON file
def getPPJSON(partid):
  with open("%s"%partid) as fp:
    data = json.load(fp)

  # need to load the refractive index data separately from the JSON
  ridata = data['ri']
  if ridata['format'] == 'gads':
    data['mList'] = [getM(path) for path in ridata['path']]
  elif ridata['format'] == 'csv':
    data['mList'] = [getMSep(path, ',') for path in ridata['path']]
  elif ridata['format'] == 'wsv':
    data['mList'] = [getMSep(path, None) for path in ridata['path']]
  else:
    print('refractive index format %s not yet supported'%ridata['format'])
    sys.exit()

  # do convenience format corrections to du psd data
  psddata = data['psd']
  if psddata['type'] == 'du':
    # see if we're using the lazy no-subbin format
    rMinMaj = []
    rMaxMaj = []
    for ai, a in enumerate(psddata['params']['rMinMaj']):
      b = psddata['params']['rMaxMaj'][ai]
      if not isinstance(a, list): # not a list, make it one
        rMinMajVal = [a]
        rMaxMajVal = [b]
      else:
        rMinMajVal = a
        rMaxMajVal = b
      
      rMinMaj.append(rMinMajVal)
      rMaxMaj.append(rMaxMajVal)

    # save new lists back
    data['psd']['params']['rMinMaj'] = rMinMaj
    data['psd']['params']['rMaxMaj'] = rMaxMaj

  return data

# Top level call, currently only works for JSON format files
def getParticleParams(partID, datatype):
  if datatype == 'json':
    return getPPJSON(partID)
  else:
    print("bad datatype: %s"%datatype)
    sys.exit()

# Specialized function to read OPAC provided table of refractive indices
def getM(fn):
  data = np.loadtxt(fn, skiprows=17, max_rows=61,comments=None,unpack=True, usecols=[1,8,9])
  data[0] *= 1e-6 # convert wavelength to meter units
  return data

# Specialized function to read CSV or WSV files of format "wavelength[um],refreal,refimag"
def getMSep(fn, sep):
  data = np.loadtxt(fn, unpack=True, delimiter=sep)
  data[0] *= 1e-6 # convert wavelength to meter units
  return data

# Specialized function to read the HITRAN provided table of water refractive indices
def getWaterM():
  data = np.loadtxt('data/refrac.water.txt', skiprows=13, unpack=True, usecols=[0,1,2])
  data[0] *= 1e-6 # convert wavelength to meter units
  return data

# Function to calculate the humidified size of particles
def humidityGrowth(params, siz0, rh, allrh): # assume all different sizes grow identically, not sure if that's true...
  rhi = list(allrh).index(rh)
  rhtype = params['type'] 
  rhp = params['params']
  # Simple growth factor specified
  if rhtype == 'simple' or rhtype == 'trivial':
    return siz0 * rhp['gf'][rhi]
  # Growth calculated after Gerber [1985]
  elif rhtype == 'ss':
     siz = siz0 * 100 # convert to cm
     if rh == 0.0:
       return siz0
     else:
       c1 = rhp['c1']
       c2 = rhp['c2']
       c3 = rhp['c3']
       c4 = rhp['c4']
       return (c1 * siz ** c2 / (c3 * siz ** c4 - np.log10(rh)) + siz ** 3.) ** (1./3.) / 100.
  elif rhtype == 'su':
    temp = rhp['temp']
    if rh == 0.0:
      return siz0
    else:
      return siz0*float(carma_utils.grow_v75(rh,siz0,temp=temp))

# Function to return a lognormal distribution in terms of size parameter and parameters
# Returns dndr assuming parameter rmode is the mode of a number distribution, or
# equivalently dNdx since it is done in size parameter space
def getLogNormPSD(rmode, sigma, xxArr, lambd, rmax, rmin):
  xconv = 2 * np.pi / lambd

  # get distribution variables in size parameter space
  xmode = rmode * xconv
  xmax = rmax * xconv
  xmin = rmin * xconv

  dNdx =   1./(xxArr * (2*np.pi) ** 0.5 * np.log(sigma) ) \
         * np.exp(-(np.log(xxArr/xmode)**2) / (2. * np.log(sigma) ** 2))

  dNdx[np.where(xxArr >= xmax)] = 0. # set values greater than xmax to zero
  dNdx[np.where(xxArr <= xmin)] = 0. # set values less than xmin to zero

  return dNdx
