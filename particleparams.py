import numpy as np
import os
import json
import sys

def getPPJSON(partid):
  with open("%s"%partid) as fp:
    data = json.load(fp)

  # need to load the refractive index data separately from the JSON
  ridata = data['ri']
  ridata2 = data['shellRi']
  
  ridatas = [ridata, ridata2]
  keys = ['mList', 'mList2']

  for iii in range(len(ridatas)):
    useridata = ridatas[iii]
    usekey = keys[iii]

    if useridata['format'] == 'gads':
      data[key] = [getM(path) for path in useridata['path']]
    elif useridata['format'] == 'csv':
      data[key] = [getMSep(path, ',') for path in useridata['path']]
    elif useridata['format'] == 'wsv':
      data[key] = [getMSep(path, None) for path in useridata['path']]
    else:
      print('refractive index format %s not yet supported'%useridata['format'])
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

def getParticleParams(partID, datatype):
  if datatype == 'json':
    return getPPJSON(partID)
  else:
    print("bad datatype: %s"%datatype)
    sys.exit()

def getM(fn):
  data = np.loadtxt(fn, skiprows=17, max_rows=61,comments=None,unpack=True, usecols=[1,8,9])
  data[0] *= 1e-6 # convert to meter units
  return data

def getMSep(fn, sep):
  data = np.loadtxt(fn, unpack=True, delimiter=sep)
  data[0] *= 1e-6 # convert to meter units
  return data

def getWaterM():
  data = np.loadtxt('data/refrac.water.txt', skiprows=13, unpack=True, usecols=[0,1,2])
  data[0] *= 1e-6 # convert to meter units
  return data

def humidityGrowth(params, siz0, rh, allrh): # assume all different sizes grow identically, not sure if that's true...
  rhi = list(allrh).index(rh)
  rhtype = params['type'] 
  rhp = params['params']
  if rhtype == 'simple' or rhtype == 'trivial':
    return siz0 * rhp['gf'][rhi]
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

def getLogNormPSD(xxArr, rmode, rmax, rmin, sigma, lambd):
  xconv = 2 * np.pi / lambd
  # get distribution variables in x space
  xmode = rmode * xconv
  xmax = rmax * xconv
  xmin = rmin * xconv

  dist3 = 1./(xxArr * (2*np.pi) ** 0.5 * np.log(sigma) ) * np.exp(-(np.log(xxArr/xmode)**2) / (2. * np.log(sigma) ** 2))

  dist3[np.where(xxArr >= xmax)] = 0. # set values greater than xmax to zero
  dist3[np.where(xxArr <= xmin)] = 0. 

  return dist3
