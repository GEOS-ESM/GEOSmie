#!/usr/bin/env python3

import numpy as np
import netCDF4 as nc
import os
from scipy.interpolate import interp1d

PFX = os.path.join('/misc', 'opk01', 'netcdf-integ')

def read_solar_flux(useWavenum):
  fn = 'solar_flux_TOA_1cm-1.txt'
  data = np.loadtxt(os.path.join('data', 'solarflux', fn))
  if not useWavenum:
    # convert from wavenumber to wavelength 
    data[:,0] = 1. / data[:,0] * 0.01
  return data

def doAverage(lam, varIn, bandl, bandr, useWavenum, solardata):
  # TODO! solardata is not used in any way even if the parameter is set to true
  if solardata is None:
    useSolar = False
  else:
    useSolar = True

  if not useWavenum: 
    # get wavenumbers from wavelengths in bandlow/up
    bandl = (100. * bandl) ** (-1)
    bandr = (100. * bandr) ** (-1)

  wnum = (100.*lam)**(-1)  # wavelength [m] to wavenumber [cm]
  # wavenumbers of bands are in integer units of [cm-1]
  num_subbin = 100
  if bandl < bandr:
    beg = bandl
    end = bandr
  else:
    beg = bandr
    end = bandl

  wnum_interpolate = np.linspace(beg, end, num_subbin)

  # ensure we don't exceed the available wavenumbers
  wnum_interpolate = np.clip(wnum_interpolate, np.min(wnum), np.max(wnum))
  
  #  Interpolate the input to the wnum_interpolate
  f_ipol = interp1d(wnum,varIn)
  varIn_interpolate = f_ipol(wnum_interpolate)
  #  Linear weighting yields...
  varOut = np.sum(varIn_interpolate) / len(varIn_interpolate)
  return varOut

def doSolarWeighting(varIn, bandl, bandr, solar):
  num_subbin = 100
  if bandl < bandr:
    beg = bandl
    end = bandr
  else:
    beg = bandr
    end = bandl

  wnum_interpolate = np.linspace(beg, end, num_subbin)
  solar_wnum = solar[:,0]
  solar_flux = solar[:,1]
  
  # ensure we don't exceed the available wavenumbers
  wnum_interpolate = np.clip(wnum_interpolate, np.min(solar_wnum), np.max(solar_wnum))

  f_ipol = interp1d(solar_wnum, solar_flux)
  wght = f_ipol(wnum_interpolate)

  varOut = np.sum(varIn_interpolate*wght) / np.sum(wght)

def getBands(mode):
  useWavenum = True # use wavenumbers by default
  numBandsMod = 0 # by default band numbers are the same as assigned here
  if mode == 'GEOS5':
    useWavenum = False
    lBandLow = np.array([.175, .225, .285, .300, .325, .400, .690, 1.220, 2.270,\
                29.412, 18.519, 12.5, 10.204, 9.091, 8.230, 7.246, 5.263, 3.333, 16.129])*1.e-6
    lBandUp  = np.array([.225, .285, .300, .325, .400, .690, 1.220, 2.270, 3.850,\
                40., 29.412, 18.519, 12.5, 10.204, 9.091, 8.230, 7.246, 5.263, 18.519])*1.e-6
    numBandsMod = -1 # reduce the output dimension by one due to special combination of two bands

  elif mode == 'RRTMG':
    wvnl_sw = [2600., 3250., 4000., 4650., 5150., 6150., 7700., 8050.,\
               12850., 16000., 22650., 29000., 38000., 820.]
    wvnr_sw = [3250., 4000., 4650., 5150., 6150., 7700., 8050.,\
               12850., 16000., 22650., 29000., 38000., 50000., 2600.]
    wvnl_lw = [10., 350., 500., 630., 700., 820., 980., 1080.,\
               1180., 1390., 1480., 1800., 2080., 2250., 2380., 2600.]
    wvnr_lw = [350., 500., 630., 700., 820., 980., 1080., 1180.,\
               1390., 1480., 1800., 2080., 2250., 2380., 2600., 3250.]
    lBandLow = wvnl_sw + wvnl_lw
    lBandUp  = wvnr_sw + wvnr_lw

  elif mode == 'RRTMGP':
    wvnl_sw = [820.,2680., 3250., 4000., 4650., 5150., 6150., 7700., 8050., \
               12850., 16000., 22650., 29000., 38000.]
    wvnr_sw = [2680.,3250., 4000., 4650., 5150., 6150., 7700., 8050., \
               12850., 16000., 22650., 29000., 38000., 50000.]
    # LW - longwave bands
    wvnl_lw = [10., 250., 500., 630., 700., 820., 980., 1080., \
               1180., 1390., 1480., 1800., 2080., 2250., 2390., 2680.]
    wvnr_lw = [250., 500., 630., 700., 820., 980., 1080., 1180., \
             1390., 1480., 1800., 2080., 2250., 2390., 2680., 3250.]
    lBandLow = wvnl_sw + wvnl_lw
    lBandUp  = wvnr_sw + wvnr_lw

  elif mode == 'PURDUE':
    useWavenum = False
    lBandLow = np.array([.175, .225, .245, .280, .295, .310, .325, .400, .700,\
                1.220, 2.270, 3.33, 5.26, 7.25, 8.23, 9.09, 10.2, 12.5,\
                16.13, 18.52, 29.41])*1.e-6
    lBandUp  = np.array([.225, .280, .260, .295, .310, .320, .400, .700,\
                1.220, 2.270, 10.00, 5.26, 7.25, 8.23, 9.09, 10.2,\
                12.5, 18.52, 18.52, 29.41, 40.0])*1.e-6

# Compute the mean value of the band range in [m]
  bandMeanM = [(x+y)/2. for x, y in zip(lBandLow,lBandUp)]
  if useWavenum:
    bandMeanM = [(x * 100.)**(-1) for x in bandMeanM]

  numBands = len(lBandLow) + numBandsMod
  return lBandLow, lBandUp, bandMeanM, useWavenum, numBands

def fun(data, part, opfn, mode, useSolar, noIR):
  solardata = None
  if useSolar:
    solardata = read_solar_flux(True)
  #print(data.variables.keys())

  try:
    lam = data.variables['wavelength'][:]
    oppclassic = False
    radiusNm = 'bin'
    lamNm = 'wavelength'
  except:
    lamNm = 'lambda'
    lam = data.variables[lamNm][:]
    oppclassic = True
    radiusNm = 'radius'
    print("Operating on a legacy file")
  rh = data.variables['rh'][:]
  reff = data.variables['rEff'][:]
  radius = data.variables[radiusNm][:]

  lBandLow, lBandUp, bandMeanM, useWavenum, nbands = getBands(mode)

  nbands_orig = len(lBandLow) # this is the dimension we do calculations at
  nrh = len(rh)
  nbin = len(reff[:,0])

  output = {}

  varsToAverage = ['qsca', 'qext', 'bsca', 'bext', 'g', 'bbck', 'refreal', 'refimag']

  for varName in varsToAverage:
    if oppclassic:
      output[varName] = np.zeros([nbin, nrh, nbands_orig])
    else:
      output[varName] = np.zeros([nbin, nbands_orig, nrh])

  for irh in range(0, nrh):
    for ilam in range(0, nbands_orig):
      for ibin in range(0, nbin):
        # average each key quantity here...
        for varName in varsToAverage:
          if oppclassic:
            varIn = data.variables[varName][ibin,irh,:]
          else:
            varIn = data.variables[varName][ibin,:,irh]
          ret = doAverage(lam, varIn, lBandLow[ilam], lBandUp[ilam], useWavenum, useSolar)

          # save to output
          if oppclassic:
            output[varName][ibin, irh, ilam] = ret
          else:
            output[varName][ibin, ilam, irh] = ret


  # if noIR then set IR band values to zero/epsilon
  if(noIR and mode == 'GEOS5'):
    ind = np.where(lBandLow > 3e-6)
    if oppclassic:
      output['qsca'][:, :, ind] = 0
      output['bsca'][:, :, ind] = 0
      output['qext'][:, :, ind] = 1e-32
      output['bext'][:, :, ind] = 1e-32
    else:
      output['qsca'][:, ind, :] = 0
      output['bsca'][:, ind, :] = 0
      output['qext'][:, ind, :] = 1e-32
      output['bext'][:, ind, :] = 1e-32

  # TODO! 
  # if mode == tGEOS5 we need to average bands 0 and 2 together and set that as
  # the lowest band, then move everything down by one?
  # above is probably not correct...I'm not sure what we're actually supposed to do with this
  # I think the idea is that we average band 0 with band 2, then just kick band 0 out
  # so then the bands will be [band1, (band0+band2)/2, band3, band4, ...)
  # I have no clue why this is done...better check with Pete
  if mode == 'GEOS5':
    for key, val in output.items():
      # this weighing needs to be done in wavelength space
      bw1 = lBandUp[0] - lBandLow[0]
      bw2 = lBandUp[2] - lBandLow[2]
      if oppclassic:
        avg = (bw1 * val[:, :, 0] + bw2 * val[:, :, 2]) / (bw1 + bw2)
        # assign the average to index 2
        val[:, :, 2] = avg
        # drop the index 0 and assign back to the ret dict
        output[key] = val[:, :, 1:]
      else:
        avg = (bw1 * val[:, 0, :] + bw2 * val[:, 2, :]) / (bw1 + bw2)
        # assign the average to index 2
        val[:, 2, :] = avg
        # drop the index 0 and assign back to the ret dict
        output[key] = val[:, 1:, :]

      bandMeanM = bandMeanM[1:]

  opncdf = nc.Dataset(opfn, 'w')
  # create dimensions
  opncdf.createDimension('rh',nrh)
  opncdf.createDimension(lamNm,nbands)
  opncdf.createDimension(radiusNm,nbin)
  nchar = 80
  opncdf.createDimension('nchar',nchar)
  for varName in output.keys():
    opncdf.createVariable(varName, 'f8', (radiusNm, lamNm, 'rh'), compression='zlib')
    #print(opncdf.variables[varName].shape)
    #print(output[varName].shape)
    #print(opncdf.variables[varName][:].shape)
    #print(output[varName].shape)
    opncdf.variables[varName][:] = output[varName]

  # special variables: lambda, rh, qname
  opncdf.createVariable(lamNm, 'f8', (lamNm))
  opncdf.variables[lamNm][:] = bandMeanM

  opncdf.createVariable('rh', 'f8', ('rh'))
  opncdf.variables['rh'][:] = data.variables['rh'][:]

  opncdf.createVariable('qname', 'c', (radiusNm, 'nchar'))
  qname = ['%s%03d'%(part, i) for i in range(1, len(radius)+1)]

  # add low and high limit information for bands
  # convert to wavelength if needed
  if mode == 'RRTMG':
    # convert from wavenumber to wavelength
    lBandLow = np.array(lBandLow) ** (-1) * 0.01
    lBandUp = np.array(lBandUp) ** (-1) * 0.01
    lBandLow0 = lBandUp[::-1]
    lBandUp0 = lBandLow[::-1]
    lBandUp = lBandUp0
    lBandLow = lBandLow0
  lBandLow = lBandLow * 1e6
  lBandUp = lBandUp * 1e6
  opncdf.createVariable('bandLow', 'f8', (lamNm))
  opncdf.variables['bandLow'][:] = lBandLow
  opncdf.variables['bandLow'].long_name = 'Lower edges of the bands'
  opncdf.variables['bandLow'].units = 'micrometers'

  opncdf.createVariable('bandUp', 'f8', (lamNm))
  opncdf.variables['bandUp'][:] = lBandLow
  opncdf.variables['bandUp'].long_name = 'Upper edges of the bands'
  opncdf.variables['bandUp'].units = 'micrometers'

  # backwards way of writing the strings
  for qi, qq in enumerate(qname):
    for ci, char in enumerate(qq):
      opncdf.variables['qname'][qi, ci] = char

  radvars = [radiusNm]
  for var in radvars:
    opncdf.createVariable(var, 'f8', (radiusNm))
    opncdf.variables[var][:] = data.variables[var][:]

  radvars2 = ['rLow', 'rUp']
  # need to drop rh dimension for these variables
  for var in radvars2:
    opncdf.createVariable(var, 'f8', (radiusNm))
    opncdf.variables[var][:] = data.variables[var][:, 0]

  # variables that have dimensions radius, rh (rEff, rMass)
  radrhvars = ['rEff', 'rMass']
  for var in radrhvars:
    opncdf.createVariable(var, 'f8', (radiusNm, 'rh'))
    opncdf.variables[var][:] = data.variables[var][:]

  opncdf.close()

    
def printstuff():
  modes = ['GEOS5', 'RRTMG']
  printbands = []
  for mode in modes:
    lBandLow, lBandUp, useWavenum, nbands = getBands(mode)
    if mode == 'RRTMG':
      # convert from wavenumber to wavelength 
      #data[:,0] = 1. / data[:,0] * 0.01
      lBandLow = np.array(lBandLow) ** (-1) * 0.01
      lBandUp = np.array(lBandUp) ** (-1) * 0.01
      lBandLow0 = lBandUp[::-1]
      lBandUp0 = lBandLow[::-1]
      lBandUp = lBandUp0
      lBandLow = lBandLow0
    print(mode)
    lBandLow = lBandLow * 1e6
    lBandUp = lBandUp * 1e6
    for i in range(len(lBandLow)):
      print('%.4f - %.4f'%(lBandLow[i], lBandUp[i]))
    print()

def main():
  #printstuff()
  partnames = ['bc', 'oc', 'su', 'ss', 'du', 'ni', 'brc', 'du', 'du-lognorm']
  #partnames = ['bc']#, 'oc', 'su', 'ss', 'du', 'ni', 'brc']
  partnames = ['du-lognorm']
  for part in partnames:
    print('Start particle %s'%part)
    processParticle(part)

def processFileForBandMode(filepath, partcode, opdir, bandmode, useSolar, noIR):
  data = nc.Dataset(filepath, 'r')

  # get base file name
  fn = os.path.basename(filepath)

  # replace .nc with .BANDMODE.nc
  opfn = fn.replace('.nc', '.%s.nc'%bandmode)
  oppath = os.path.join(opdir, opfn)

  fun(data, partcode, oppath, bandmode, useSolar, noIR)

def processParticle(part):
  # separate function for each particle so we don't have to read in the netCDF file
  # for every iteration of the band
  fn = 'integ-%s-raw.nc'%part
  fullfn = os.path.join(PFX, fn)

  useSolar = False
  noIR = False
  solarpart = ''
  if useSolar:
    solarpart = '.solar'
  irpart = ''
  if noIR:
    irpart = '.noir'
  fnstrip = 'bands-%s'%(part)
  pfx = os.path.join('/misc', 'opk01', 'geos-bands')

  modes = ['GEOS5', 'RRTMG', 'RRTMGP', 'PURDUE']
  for mode in modes:
    print('STARTING MODE %s'%mode)
    opfn = '%s.%s%s%s.nc'%(fnstrip, mode, solarpart, irpart)
    oppath = os.path.join(pfx, opfn)
    processFileForBandMode(fullfn, part, oppath, mode, useSolar, noIR)

if __name__ == '__main__':
  main()
