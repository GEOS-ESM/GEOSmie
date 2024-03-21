"""
Read in a NetCDF file containing legendre-expanded phase matrix elements
and write them out in a physical format at angles of specific gaussian nodes,
then feed this table into Mishchenko's spher_expand.f code that calculates
the VLIDORT-compatible expansion

NOTE! Pete's legendre coefficients should be evaluated with homebrew legendre code,
for some reason the standard Python legendre code doesn't work (possibly different
normalization etc)

General working concept:
-open a netcdf file
-loop over all lambda
-loop over all RH
-read in s11, s12, etc, write them in a temp table file, run Mishchenko's code for that file,
read in Mishchenko's output, save that output back to NetCDF under another variable name
(or separate file)
-delete temporary file and Mischenko's output
"""

import netCDF4
import os
import sys
import numpy as np
from numpy.polynomial.legendre import leggauss
import numba
import multiprocessing
import time
import shutil

# x is cos(ang)
@numba.jit(nopython=True)
def hbleg(x, pmoms, leg):
  leg[0,:] = 1.
  leg[1,:] = x
  for imom in range(2,len(pmoms)):
    leg[imom] =  1. / imom*((2.*imom-1.)*x*leg[imom-1]-(imom-1.)*leg[imom-2])

  # better to do with numpy dot
  ret = np.zeros_like(x)
  for nm in range(len(pmoms)):
    ret += pmoms[nm] * leg[nm]
  return ret


def createVariablesLegendre(ncdf, numExpand, numangs, linangs):
  ncdf.createVariable('pmom2', 'f8', ('nPol', 'nMom', 'radius', 'rh', 'lambda')) # we copy the original moments here
  ncdf.createDimension('nMom-GSF', numExpand)

  ncdf.variables['pmom2'].long_name = ncdf.variables['pmom'].long_name
#  ncdf.variables['pmom'].long_name =\
#'nts of the generalized spherical functions, ordered over nGSFun as Alpha1, Alpha2, Alpha3, Alpha4, Beta1, Beta2 (Mishchenko)\
#ota, Alpha, Zeta, Delta, Gamma, Epsilon (VLIDORT)'
  #del ncdf.variables['pmom']

  #ncdf.createVariable('pmom', 'f8', ('nPol', 'nMom-GSF', 'radius', 'rh', 'lambda')) # we copy the original moments here
  ncdf.variables['pmom'].long_name =\
'Mts of the generalized spherical functions order over nPol as the 11, 12, 33, 34, 22, 44 \
elts of the scattering law matrix for a Fourier decomposition of the phase matrix such as those used in discrete ordinate solutions to the RTE'

  ncdf.createDimension('scattering_angle', numangs)
  ncdf.createVariable('scattering_angle', 'f8', ('scattering_angle'))
  ncdf.variables['scattering_angle'][:] = linangs
  ncdf.variables['scattering_angle'].long_name = 'Scattering angle in degrees'

  ncdf.createVariable('phase_matrix', 'f8', ('nPol', 'scattering_angle', 'radius', 'rh', 'lambda')) # we copy the original moments here
  #ncdf.variables['phase_matrix'].long_name = 'Phase matrix elements, ordered over nPol as P11, P22, P33, P44, P12, P34'
  ncdf.variables['phase_matrix'].long_name = 'Phase matrix elements, ordered over nPol as P11, P12, P33, P34, P22, P44'

def createVariablesPyGeosMie(ncdf, numExpand):
  ncdf.createDimension('m', numExpand)
  ncdf.createVariable('pmom', 'f8', ('bin','wavelength','rh','p','m'), compression='zlib') # we copy the original moments here
  ncdf.variables['pmom'].long_name =\
'Moments of the generalized spherical functions order over nPol as the 11, 12, 33, 34, 22, 44 \
elements of the scattering law matrix for a Fourier decomposition of the phase matrix such as those used in discrete ordinate solutions to the RTE'

def createVariablesIce(ncdf, numExpand, numangs, linangs):
  ncdf.createDimension('rh', 1)
  ncdf.createVariable('rh', 'f8', ('rh'))

  # rename old pmom to pmom_orig since it is missing the rh dimension
  ncdf.renameVariable('pmom', 'pmom_orig')

  # copy pmom_orig content to pmom2 for later validation
  ncdf.createVariable('pmom2', 'f8', (keydic['nPol'], keydic['nMom'], keydic['nradius'], 'rh', keydic['nlambda'])) # we copy the original moments here

  # create a new pmom variable that will hold the generalized spherical functions
  ncdf.createVariable('pmom', 'f8', (keydic['nPol'], keydic['nMom'], keydic['nradius'], 'rh', keydic['nlambda']))

  # set pmom variables from pmom_orig
  # do this in a loop since the dimensions are not aligned/correctly ordered
  for lami, lam in enumerate(ncdf.variables[keydic['lambda']]):
    for rhi, rh in enumerate(ncdf.variables['rh']):
      for radind, radius in enumerate(ncdf.variables[keydic['radius']]):
        for poli in range(ncdf.dimensions['nphamat'].size):
          ncdf.variables['pmom'][poli, :, radind, rhi, lami] = ncdf.variables['pmom_orig'][lami, radind, poli, :]

  ncdf.createDimension('nMom-GSF', numExpand)

  ncdf.createVariable('lambda', 'f8', (keydic['nlambda']))
  ncdf.createVariable('radius', 'f8', (keydic['nradius']))

  ncdf.variables['lambda'][:] = ncdf.variables[keydic['lambda']][:]
  ncdf.variables['radius'][:] = ncdf.variables[keydic['radius']][:]

  ncdf.variables['pmom2'].long_name = ncdf.variables['pmom_orig'].long_name
  ncdf.variables['pmom'].long_name =\
'Moments of the generalized spherical functions order over nPol as the 11, 12, 33, 34, 22, 44 \
elements of the scattering law matrix for a Fourier decomposition of the phase matrix such as those used in discrete ordinate solutions to the RTE'

  ncdf.createDimension('scattering_angle', numangs)
  ncdf.createVariable('scattering_angle', 'f8', ('scattering_angle'))
  ncdf.variables['scattering_angle'][:] = linangs
  ncdf.variables['scattering_angle'].long_name = 'Scattering angle in degrees'

  ncdf.createVariable('phase_matrix', 'f8', (keydic['nPol'], 'scattering_angle', keydic['nradius'], 'rh', keydic['nlambda'])) # we copy the original moments here
  ncdf.variables['phase_matrix'].long_name = 'Phase matrix elements, ordered over nPol as P11, P12, P33, P34, P22, P44'

  ncdf.createVariable('pback', 'f8', (keydic['nradius'], 'rh', keydic['nlambda'], keydic['nPol']))
  ncdf.createVariable('bbck', 'f8', (keydic['nradius'], 'rh', keydic['nlambda']))
  ncdf.createVariable('g', 'f8', (keydic['nradius'], 'rh', keydic['nlambda']))
  ncdf.createVariable('bsca', 'f8', (keydic['nradius'], 'rh', keydic['nlambda']))
  ncdf.createVariable('bext', 'f8', (keydic['nradius'], 'rh', keydic['nlambda']))
  ncdf.createVariable('qsca', 'f8', (keydic['nradius'], 'rh', keydic['nlambda']))
  ncdf.createVariable('qext', 'f8', (keydic['nradius'], 'rh', keydic['nlambda']))

  ncdf.variables['g'].long_name = 'Asymmetry parameter (UNVALIDATED)'
  ncdf.variables['pback'].long_name = 'Backscattering phase function (UNVALIDATED)'
  ncdf.variables['bbck'].long_name = 'Backscattering mass efficiency (UNVALIDATED)'
  ncdf.variables['bsca'].long_name = 'Scattering mass efficiency (UNVALIDATED)'
  ncdf.variables['bext'].long_name = 'Extinction mass efficiency (UNVALIDATED)'
  ncdf.variables['qsca'].long_name = 'Scattering efficiency (UNVALIDATED)'
  ncdf.variables['qext'].long_name = 'Extinction efficiency (UNVALIDATED)'

  ncdf.variables['rh'][0] = 0.0

def createVariablesGRASP(ncdf, numExpand):
  ncdf.createDimension('nPol', 6)
  ncdf.createDimension('nMom', numExpand)
  ncdf.createDimension('radius', len(ncdf.variables['sizeBin']))
  ncdf.createDimension('rh', 1)

  ncdf.createVariable('rh', 'f8', ('rh'))
  ncdf.createVariable('radius', 'f8', ('radius'))
  ncdf.createVariable('pback', 'f8', ('radius', 'rh', keydic['lambda'], keydic['nPol']))
  ncdf.createVariable('bbck', 'f8', ('radius', 'rh', keydic['lambda']))
  #ncdf.createVariable('lidar_ratio_debug', 'f8', ('radius', 'rh', keydic['lambda']))
  ncdf.createVariable('g', 'f8', ('radius', 'rh', keydic['lambda']))
  ncdf.createVariable('bsca', 'f8', ('radius', 'rh', keydic['lambda']))
  ncdf.createVariable('bext', 'f8', ('radius', 'rh', keydic['lambda']))
  ncdf.createVariable('qsca', 'f8', ('radius', 'rh', keydic['lambda']))
  ncdf.createVariable('qext', 'f8', ('radius', 'rh', keydic['lambda']))

  ncdf.variables['radius'][:] = ncdf.variables[keydic['rv']][:] * 1e-6
  ncdf.variables['g'].long_name = 'Asymmetry parameter (UNVALIDATED)'
  ncdf.variables['pback'].long_name = 'Backscattering phase function (UNVALIDATED)'
  ncdf.variables['bbck'].long_name = 'Backscattering mass efficiency (UNVALIDATED)'
  ncdf.variables['bsca'].long_name = 'Scattering mass efficiency (UNVALIDATED)'
  ncdf.variables['bext'].long_name = 'Extinction mass efficiency (UNVALIDATED)'
  ncdf.variables['qsca'].long_name = 'Scattering efficiency (UNVALIDATED)'
  ncdf.variables['qext'].long_name = 'Extinction efficiency (UNVALIDATED)'

  ncdf.variables['rh'][0] = 0.0
  ncdf.createVariable('pmom', 'f8', ('nPol', 'nMom', 'radius', 'rh', keydic['lambda']))

def convertData(ncdf, mode, ice, whichproc, radind, rhi, lami, rhop0, num_gauss, linangs):
  linvals = np.zeros([7,len(linangs)])
  linvals[0,:] = linangs
  if mode == 'pygeos':
    ang = ncdf.variables['ang'][:]
    numang = len(ang) 
    allvals = np.zeros([7,numang])
    keys = ['s11', 's22', 's33', 's44', 's12', 's34'] # order Mischenko's code expects
    allvals[0,:] = ang
    for ki, key in enumerate(keys):
      allvals[ki+1, :] = ncdf.variables[key][radind, lami, rhi, :]

    # write the temp file
    tempfn = 'tempfile%d.txt'%whichproc
    np.savetxt(tempfn, allvals.T)
    os.system('./a.out %s > /dev/null'%tempfn)
    newdata = np.loadtxt('%s.expan_coeff'%tempfn, skiprows=1, unpack=True)
  elif mode == 'legendre':
    gpoints, gweights = leggauss(num_gauss)
    gpoints = gpoints[::-1]
    gpointsdeg = np.degrees(np.arccos(gpoints))

    pmoms = ncdf.variables['pmom']
    allvals = np.zeros([7,num_gauss])

    allvals[0,:] = gpointsdeg # code expects angle input in degrees
    # mischenko expects the order p11, p22, p33, p44, p12, p34
    #pete's order: P11, P12, P33, P34, P22, P44
    # ice files have the same order
    indexorder = [0,4,2,5,1,3]

    # we need to evaluate the legendre-based phase matrix values to angle-based for the fortran code
    for ii, npol in enumerate(indexorder):
      moms = pmoms[npol,:,radind,rhi,lami]
      leg = np.zeros([len(moms), len(gpoints)], dtype=np.float64)
      vals = hbleg(gpoints, moms, leg)
      allvals[1+ii,:] = vals

      leg2 = np.zeros([len(moms), len(linangs)], dtype=np.float64)
      vals2 = hbleg(np.cos(np.radians(linangs)), moms, leg2)
      linvals[1+ii,:] = vals2 # mischenko's order

    # write the temp file
    tempfn = 'tempfile%d.txt'%whichproc
    np.savetxt(tempfn, allvals.T)
    os.system('./a.out %s > /dev/null'%tempfn)
    newdata = np.loadtxt('%s.expan_coeff'%tempfn, skiprows=1, unpack=True)

    if ice: # for ice we have to calculate extra variables
      """
      Calculate pback
      """
      p11 = ncdf.variables['phase'][lami, radind, 0, :]
      for ii, npol in enumerate([0,4,2,5,1,3]):
        ncdf.variables['pback'][radind, rhi, lami, npol] = allvals[ii+1][-1] # this should be the backscattering angle
      ncdf.variables['bext'][radind,rhi,lami] = ncdf.variables['ext'][lami,radind]
      ncdf.variables['bsca'][radind,rhi,lami] = ncdf.variables['ext'][lami,radind] * ncdf.variables['ssa'][lami, radind]
      pBck = p11[-1] # p11 at backscatter direction
      qBck = pBck * ncdf.variables['qsca'][radind,rhi,lami] / (4*np.pi)
      rho = ncdf.variables['rho'][0]
      reff = ncdf.variables[keydic['radius']][radind]
      bBck = 3./4.*qBck / rho / reff
      ncdf.variables['bbck'][radind,rhi,lami] = bBck * 1e6

      """
      Calculate g
      """

      angs = ncdf.variables['theta'][lami, radind, 0, :]
      angs = np.radians(angs)
      p11n = p11 / np.sum(p11 * np.sin(angs))
      g = np.sum(np.cos(angs) * np.sin(angs) * p11n)

      ncdf.variables['g'][radind, rhi, lami] = g

  else: # physical space, i.e. GRASP
    numang = 181 # hardcoded in GRASP
    allvals = np.zeros([7,numang])
    angs = np.linspace(0.,180.,numang)
    allvals[0,:] = angs
    keys = ['p11', 'p22', 'p33', 'p44', 'p12', 'p34'] # order Mischenko's code expects
    for ki, key in enumerate(keys):
      allvals[ki+1, :] = ncdf.variables[key][radind, lami, :]

    p11 = allvals[1][:] # p11 at backscatter direction

    """
    Calculate pback
    """
    for ii, npol in enumerate([0,4,2,5,1,3]):
      ncdf.variables['pback'][radind, lami, rhi, npol] = allvals[ii+1][-1] # this should be the backscattering angle

    """
    Calculate bbck etc
    """

    convfact = 1. / (rhop0[radind] * 1e-6)
    foo1 = ncdf.variables['bsca_vol'][radind,:]
    ncdf.variables['bsca'][radind,lami,rhi] = ncdf.variables['bsca_vol'][radind,lami] * convfact
    ncdf.variables['bext'][radind,lami,rhi] = ncdf.variables['bext_vol'][radind,lami] * convfact
    ncdf.variables['qsca'][radind,lami,rhi] = ncdf.variables['bsca_vol'][radind,lami] * 4./3. * ncdf.variables['rEff'][radind]
    ncdf.variables['qext'][radind,lami,rhi] = ncdf.variables['bext_vol'][radind,lami] * 4./3. * ncdf.variables['rEff'][radind]
    pBck = ncdf.variables['p11'][radind,lami,-1] # p11 at backscatter direction
    qBck = pBck * ncdf.variables['qsca'][radind,lami,rhi] / (4*np.pi)
    bBck = 3./4.*qBck / rhop0[radind] / ncdf.variables['rEff'][radind]
    ncdf.variables['bbck'][radind,lami,rhi] = bBck * 1e6

    """
    Calculate g
    """

    angs = np.radians(angs)
    # TODO: check the definition, should be between 1 and -1
    p11n = p11 / np.sum(p11 * np.sin(angs))
    g = np.sum(np.cos(angs) * np.sin(angs) * p11n)

    ncdf.variables['g'][radind, lami, rhi] = g

    # write the temp file
    tempfn = 'tempfile%d.txt'%whichproc
    np.savetxt(tempfn, allvals.T)
    os.system('./a.out %s > /dev/null'%tempfn)
    newdata = np.loadtxt('%s.expan_coeff'%tempfn, skiprows=1, unpack=True)

  return newdata, linvals

def fun(fn, whichproc, rhop0, ice, keydic):
  print("Processing file %s"%fn)

  infile = fn

  outdir = os.path.join('.')

  #varname = 'spher_expand_coeff' # eventual(?) name
  varname = 'pmom' # re-use the same name to let VLIDORT re-use code

  if 'GRASP' not in fn:
    if 'integ' in fn:
      mode = 'pygeos'
    else:
      mode = 'legendre'
  else:
    mode = 'physical'

  processFileRaw(infile, outdir, whichproc, rhop0, mode, ice)

def processFileRaw(infile, outdir, whichproc, rhop0, mode, ice):
  # mischenko expects the order p11, p22, p33, p44, p12, p34
  #numExpand = 2001 # this needs to match what is set in params.h
  #numExpand = 300 
  numExpand = 129 
  num_gauss = 960 # arbitrary? perhaps try adaptive
  # mishchenko's code gives errors in num_gauss is bigger than 1000, though this can probably be changed in params.h

  fn = os.path.basename(infile)
  sfx = fn.split('.')[-1] # file suffix
  outfn = fn.replace(sfx, 'GSFun.%s'%sfx)
  if numExpand != 2001:
    outfn = outfn.replace('.GSFun.%s'%sfx, '.GSFun-%d.%s'%(numExpand, sfx))

  outfile = os.path.join(outdir, outfn)

  shutil.copyfile(infile, outfile) # copy the original to the to-be-modified file
  ncdf = netCDF4.Dataset(outfile, 'r+')

  
  numangs = 181 # for linearly evaluated phase matrix
  linangs = np.linspace(0,180,numangs)

  # start creating variables required for the tables
  print('mode %s'%mode)
  if (mode == 'legendre') and not ice:
    createVariablesLegendre(ncdf, numExpand, numangs, linangs)
    keydic = {'nradius': 'radius', 'nPol': 'nPol', 'nlambda': 'lambda', 'nMom': 'nMom', 'lambda': 'lambda', 'radius': 'radius'}
  elif mode =='legendre' and ice:
    createVariablesIce(ncdf, numExpand, numangs, linangs)
    keydic = {'nradius': 'nreff', 'nPol': 'nphamat', 'nlambda': 'nlam', 'nMom': 'nmommax', 'lambda': 'wavelen', 'radius': 'reff'}
  elif (mode == 'pygeos'):
    createVariablesPyGeosMie(ncdf, numExpand)
    keydic = {'wavelength': 'wavelength', 'bin': 'bin'}
  else: # grasp
    createVariablesGRASP(ncdf, numExpand)
    keydic = {'rv': 'rv', 'nPol': 'nPol', 'lambda': 'lambda', 'radius': 'sizeBin'}

  alllambda = ncdf.variables[keydic['wavelength']]

  print(ncdf.variables.keys())
  for lami, lam in enumerate(alllambda):
    if lami % 10 == 0:
      print("lami %d of %d"%(lami+1, len(alllambda)))
    for rhi, rh in enumerate(ncdf.variables['rh']):
      for radind, radius in enumerate(ncdf.variables[keydic['bin']]):
        newdata, linvals = convertData(ncdf, mode, ice, whichproc, radind, rhi, lami, rhop0, num_gauss, linangs)
            
        """
        Start writing output variables to the file
        """
        for ii, npol in enumerate([0,4,2,5,1,3]):
          if mode == 'legendre':
            ncdf.variables['pmom2'][radind, lami, rhi, npol, :] = ncdf.variables['pmom'][radind, lami, rhi, npol, :]
          pmomdata = newdata[1+ii]
          lenpmom = len(ncdf.variables['pmom'][radind, lami, rhi, npol, :])
          pmomdata2 = np.zeros(lenpmom)
          pmomdata2[:len(pmomdata)] = pmomdata # zero-padded array
            
          ncdf.variables['pmom'][radind, lami, rhi, npol, :] = pmomdata2
          if mode == 'legendre':
            ncdf.variables['phase_matrix'][radind, lami, rhi, npol, :] = linvals[1+ii]
            angs = ncdf.variables['scattering_angle']

    if lami + 1 == len(alllambda):
      print("%s done"%fn)

def main():
  numprocmax = 1

  allfn = [\
  'optics_BC.v5_7.nc',\
  'optics_SU.v5_7.nc',\
  'optics_BRC.v12_7.nc',\
  'optics_OC.v12_7.nc',]
 # 'optics_SS.v3_7.nc']
  allfn = ['optics_SS.v3_7.nc']
  allfn = ['ic.column_8elements.050.1.cdf', 'ic.column_8elements.050.2.cdf']
  allfn = ['ic.solid_bullet_rosette.000.1.cdf','ic.solid_bullet_rosette.000.1.cdf',\
  'ic.solid_bullet_rosette.003.1.cdf', 'ic.solid_bullet_rosette.003.2.cdf',\
  'ic.solid_bullet_rosette.050.1.cdf', 'ic.solid_bullet_rosette.050.2.cdf']

  allfn = [\
  "integ-du_sphere-raw.nc",\
  "integ-bc-raw.nc",\
  "integ-brc-raw.nc",\
  "integ-oc-raw.nc",\
  "integ-su-raw.nc",\
  "integ-ni-raw.nc",\
  "integ-ss-raw.nc",\
  ]

  allfn = [\
  "integ-ni-raw.nc",\
  "integ-ss-raw.nc",\
  ]

  #allfn = ['GRASP_LUT-DUST_V4.nc']
  #allfn = ['GRASP_LUT-DUST_V3.nc']
  ##allfn = ['GRASP_LUT-DrySU_V1.nc']
  #allfn = ['GRASP_LUT-DUST_V4.nc', 'GRASP_LUT-DrySU_V1.nc']
  #allfn = ['GRASP_LUT-DrySU_V1.nc']
  #allfn = ['GRASP_LUT-DUST_V4.nc']
  #allfn = ['GRASP_LUT-DUST_V5_80nmTO20000nm.nc', 'GRASP_LUT-DUST_V5_80nmTO500nm.nc']
  #allfn = ['optics_DU.v15_6.nc']
  #allfn = ['optics_DU.v18_6.nc']
  dustrho = np.array([2500., 2650., 2650., 2650., 2650.])
  surho = np.array([1700.])
  icerho = np.array([1000.])
  

  allrho = [dustrho] # make this match the filenames
  #allrho = [surho] # make this match the filenames
  #allrho = [icerho] 

  ice = False
  #ice = False

  numproc = min(len(allfn), numprocmax)
  beg = time.time()
  if numproc > 1:
    mypool = multiprocessing.Pool(processes=numproc)
    pool_results = np.empty([len(allfn)], dtype=multiprocessing.pool.ApplyResult)

    for fni, fn in enumerate(allfn):
      thisrho = None
      if 'DUST' in fn:
        thisrho = dustrho
      elif 'DrySU' in fn:
        thisrho = surho
      this_result = mypool.apply_async(fun, (fn, fni, thisrho))
      pool_results[fni] = this_result

    for fni, fn in enumerate(allfn):
      ret = pool_results[fni].get()
  else:
    for fni, fn in enumerate(allfn):
      thisrho = None
      if 'DUST' in fn:
        thisrho = dustrho
        print('using dust rho')
      elif 'DrySU' in fn:
        thisrho = surho
        print('using su rho')
      elif 'ic' in fn:
        thisrho = icerho    
        print('using ice rho')
      else:
        print('rho not defined')
        thisrho = None
        #sys.exit()
      fun(fn, fni, thisrho, ice, {})

  end = time.time()
  print("took %.2f seconds on %d processors"%(end-beg, numproc))

def convertFile(filepath, opdir, mode, rhop0):
  # get base file name
  fn = os.path.basename(filepath)

  if mode in ['pygeos', 'legendre']:
    rhop0 = None # unset this to avoid any accidental use

  ice = False 
  # ice is handled as a subcategory of legendre so we set a separate variables
  if mode == 'ice':
    mode = 'legendre'
    ice = True

  whichproc = 0 # parallelization is not implemented for rungsf interface
  processFileRaw(filepath, opdir, whichproc, rhop0, mode, ice)

if __name__ == '__main__':
  main()
  #fun('optics_BC.v1_6.nc')
  #fun('optics_BC.v2_6.nc')
  
