#!/usr/bin/env python3
#

"""

   Implements API to access Mie LUTs for a single specie.

"""

import numpy as np
import xarray as xr
import pandas as pd

__VERSION__ = '0.9.0'

valid_names = ['aot',          'ssa',     'gf',      'gasym',   'g',   'growth_factor',
               'RefIndex',     'pmom',    'area',    'volume',  'p11', 'p22', 'pback',
               'rhod',         'rhop',    'rEff',    'bbck',    'tau',
               'bsca',         'bext',    'refreal', 'refimag',
               'aot_ssa_pmom',
               'aot_ssa_gasym' ] 

alias_names = {'gf'   : 'growth_factor',
               'tau'  : 'aot',
               'gasym': 'g'}


class MieTABLE(object):

   def __init__ ( self, filename):
      """
      Uses xarray to load Mie table from file name. On inout
      filename    --- Mie Table file name
      # Single Scattering properties (fix order of dimension, what is below comes from Fortran
      # ----------------------------
      #wavelengths     # (c) wavelengths        [m]
      #rh              # (r) RH values   [fraction]
      #reff            # (b,r) effective radius [m]
      #bext            # (b,r,c) bext values [m2 kg-1]
      #bsca            # (b,r,c) bsca values [m2 kg-1]
      #bbck            # (b,r,c) bbck values [m2 kg-1]
      #g               # (b,r,c) asymmetry parameter
      #p11             # (b,r,c) Backscatter phase function, index 1 
      #p22             # (b,r,c) Backscatter phase function, index 5
      #pmom            # (p,m,b,r,c) moments of phase function
      #pback           # (p,b,r,c)   moments of backscatter phase function
      #gf              # (b,r) hygroscopic growth factor
      #rhop            # (b,r) wet particle density [kg m-3]
      #rhod            # (b,r) dry particle density [kg m-3]
      #vol             # (b,r) wet particle volume [m3 kg-1]
      #area            # (b,r) wet particle cross section [m2 kg-1]
      #refr            # (b,r,c) real part of refractive index
      #refi            # (b,r,c) imaginary part of refractive index
      """ 
      self.filename    = filename
      self.mieDS       = xr.open_dataset(filename)
      self.vars_list   = self.mieDS.keys()
    
   def __getTable__(self, name, bin, wavelength=None):
      """
        get table directly from file
      """
      assert name in self.vars_list, name + ' is not found in the table '+ self.filename

      bin_ = bin - 1
      if 'lambda' in self.mieDS[name].dims:
         assert wavelength is not None, 'wavelength should be provided to get variable ' + name
         var = self.mieDS[name].isel({'radius':[bin_]}).interp({'lambda': [wavelength]})
      else:
         var = self.mieDS[name].isel({'radius':[bin_]})

      return var
#--
   def getVariable(self, name, bin, rh, wavelength=None, q_mass=None):
      """
       get Variables by the name in the list valid_names
      """
      assert name in valid_names, "getVariable does not support " + name
      if name in alias_names : name = alias_names[name]
      
      if name in self.vars_list :
         var = self.__getTable__(name, bin, wavelength=wavelength)
         var = var.interp(rh=rh)

      if name == 'aot' :
         assert q_mass is not None, 'aot needs q_mass as input'
         bext  = self.__getTable__('bext', bin, wavelength=wavelength)
         bext_ = bext.interp(rh=rh)
         var   = (bext_*q_mass).rename('aot')

      if name == "ssa":
         bext = self.__getTable__('bext', bin, wavelength=wavelength)           
         bsca = self.__getTable__('bsca', bin, wavelength=wavelength)           
         ssa  = bsca/bext
         var  = ssa.interp(rh=rh).rename('ssa')

      if name == 'volume':
         rhod = self.__getTable__('rhod', bin)
         gf   = self.__getTable__('gf',   bin)
         vol  = gf**3/rhod
         var = vol.interp(rh=rh).rename('volume')

      if name == 'area':
         rhod  = self.__getTable__('rhod', bin)
         gf    = self.__getTable__('gf',   bin)
         reff  = self.__getTable__('rEff', bin)
         vol   = gf**3/rhod
         area  = vol/(4./3.*reff)
         var   = area.interp(rh=rh).rename('area')

      if name == 'RefIndex':
         refr = self.__getTable__('refreal', bin, wavelength=wavelength)
         refi = self.__getTable__('refimag', bin, wavelength=wavelength)
         var = (refr.interp(rh=rh), refi.interp(rh=rh))

      if name == 'aot_ssa_pmom' or name == 'aot_ssa_gasym':
         assert q_mass is not None, name + 'needs q_mass as input'
         bext = self.__getTable__('bext', bin, wavelength=wavelength)
         bsca = self.__getTable__('bsca', bin, wavelength=wavelength)
         ssa  = (bsca/bext).interp(rh=rh).rename('ssa')
         aot  = (bext.interp(rh=rh) * q_mass).rename('aot')
         if 'pmom' in name:
            pmom = self.__getTable__('pmom', bin, wavelength=wavelength)
            pmom = pmom.interp(rh=rh)
            var  = (aot, ssa, pmom)
         if 'gasym' in name:
            gasym = self.__getTable__('g', bin, wavelength=wavelength)
            gasym = gasym.interp(rh=rh)
            var   = (aot, ssa, gasym)

      if name == 'p11':
         p11 = self.__getTable__('pback', bin, wavelength=wavelength)
         var = p11.isel({"nPol": [0]}).interp(rh=rh).rename('p11')

      if name == 'p22':
         p22 = self.__getTable__('pback', bin, wavelength=wavelength)
         var = p22.isel({"nPol": [4]}).interp(rh=rh).rename('p22')
 
      return var
       
#______________________________________________________________

if __name__ == "__main__":

   # Sample Mie Tables
   # -----------------
   dirn   = '/discover/nobackup/projects/gmao/share/dasilva/fvInput/ExtData/chemistry/AerosolOptics/v0.0.0/x/'
   Tables = [dirn + 'optics_DU.v7.nc', dirn + 'optics_OC.v2_3.nc']
   

   # Aerosol state (all species)
   # ---------------------------

   aer_Nv = '/css/gmao/geos-it/products/Y2023/M02/D05/GEOS.it.asm.aer_inst_3hr_glo_C180x180x6_v72.GEOS5294.2023-02-05T1200.V01.nc4'

   aer    = xr.open_dataset(aer_Nv).variables
  
   ######
   # Dust
   #####
   # ----
   table  = Tables[0]
   mie    = MieTABLE(table)#,wavelengths)
   q_mass = aer['DU003'] * aer['DELP'] / 9.81
   rh     = aer['RH']

   varnames =['tau', 'aot', 'gasym', 'bext', 'bsca', 'ssa', 'bbck', 'rEff', 'RefIndex', 'aot_ssa_gasym']

   for var_name in varnames:
      print(var_name)
      var  = mie.getVariable(var_name, 3, rh, wavelength=550e-9, q_mass=q_mass)

   #####
   # OC
   #####
   # --
   table = Tables[1]
   wavelengths = [470e-9, 550e-9, 670e-9, 870e-9]
   mie = MieTABLE(table)
   q_mass = aer['BCPHILIC'] * aer['DELP'] / 9.81
   rh = aer['RH']

   for var_name in varnames:
      print(var_name)
      for wavelength in wavelengths:
         print(wavelength)
         var  = mie.getVariable(var_name, 1, rh, wavelength=wavelength, q_mass=q_mass)
