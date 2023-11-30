#!/usr/bin/env python3
#

"""

   Implements API to access Mie LUTs for a single specie.

"""

import numpy as np
import xarray as xr
import pandas as pd

__VERSION__ = '0.9.0'

valid_names = ['aot',         'ssa',     'gf',     'gasym',  'g',   'growth_factor',
               'RefIndex',    'pmom',    'area',   'volume', 'p11', 'p22', 'pback',
               'rhod',        'rhop',    'rEff',   'bbck',
               'bsca',        'bext',    'refreal','refimag',
               'aot_ssa_pmom',
               'aot_ssa_gasym' ] 

class MieTABLE(object):

   def __init__ ( self, filename, wavelengths=None, nmom=None ):
      """
      Uses xarray to load Mie table from file name. On inout

      filename    --- Mie Table file name
      wavelengths --- desired wavelengths; if omitted all wavelengths on file
      nmon        --- number of moments for phase function; if omitted no phase
                      function needed (saves memory).

      """
      self.mieDS       = xr.open_dataset(filename)

      # Single Scattering properties (fix order of dimension, what is below comes from Fortran
      # ----------------------------
      self.wavelengths = self.mieDS.coords['lambda']     # (c) wavelengths        [m]
      self.rh          = self.mieDS.coords['rh']         # (r) RH values   [fraction]
      self.reff        = None                            # (b,r) effective radius [m]
      self.bext        = None                            # (b,r,c) bext values [m2 kg-1]
      self.bsca        = None                            # (b,r,c) bsca values [m2 kg-1]
      self.bbck        = None                            # (b,r,c) bbck values [m2 kg-1]
      self.g           = None                            # (b,r,c) asymmetry parameter
      self.p11         = None                            # (b,r,c) Backscatter phase function, index 1 
      self.p22         = None                            # (b,r,c) Backscatter phase function, index 5
      self.pmom        = None                            # (p,m,b,r,c) moments of phase function
      self.pback       = None                            # (p,b,r,c)   moments of backscatter phase function
      self.gf          = None                            # (b,r) hygroscopic growth factor
      self.rhop        = None                            # (b,r) wet particle density [kg m-3]
      self.rhod        = None                            # (b,r) dry particle density [kg m-3]
      self.vol         = None                            # (b,r) wet particle volume [m3 kg-1]
      self.area        = None                            # (b,r) wet particle cross section [m2 kg-1]
      self.refr        = None                            # (b,r,c) real part of refractive index
      self.refi        = None                            # (b,r,c) imaginary part of refractive index
      
      self.vars_list   = self.mieDS.keys()

      if 'refreal' in self.vars_list:
         self.refr     = self.mieDS.refreal   
      if 'refimag' in self.vars_list:
         self.refi     = self.mieDS.refimag
      if 'pmom' in self.vars_list:
         self.pmom     = self.mieDS.pmom
      if 'bbck' in self.vars_list:
         self.bbck     = self.mieDS.bbck
      if 'g' in self.vars_list:
         self.g        = self.mieDS.g
      if 'bext' in self.vars_list:
         self.bext     = self.mieDS.bext
      if 'bsca' in self.vars_list:
         self.bsca     = self.mieDS.bsca
      if 'rhop' in self.vars_list:
         self.rhop     = self.mieDS.rhop
      if 'rhod' in self.vars_list:
         self.rhod     = self.mieDS.rhod
      if 'rEff' in self.vars_list:
         self.reff     = self.mieDS.rEff
      if 'growth_factor' in self.vars_list:
         self.gf       = self.mieDS.growth_factor

      
#--
   def getVariable(self, name, bin, rh=None, q_mass=None, wavelength=None, channel=None):
      """
       get Variables by the name in the list valid_names
      """
      
      bin_ = bin - 1
      var  = None
      if name not in valid_names:
         print("getVariable does not support " + name)

      if name == 'gasym' : name = 'g'
      if name == 'gf'    : name = 'growth_factor'
      
      if name in self.vars_list:
         if 'lambda' in self.mieDS[name].dims:
            channel_ = channel
            if channel_ is None:
               assert wavelength is not None, "Either wavelength or channel should be provided as input"
               channel_= self.getChannel(wavelength)
               assert channel_ is not None, "Cannot find the wavelength: " + str(wavelength)

            var = self.mieDS[name].isel({'radius':[bin_], 'lambda':[channel_]})              
         else:
            var = self.mieDS[name].isel({'radius':[bin_]})
         if (rh is not None) :
            var = var.interp(rh = rh)

      if name == 'aot' :
         if q_mass is not None:
            bext =  self.getVariable('bext', bin, rh = rh, wavelength = wavelength, channel= channel)
            var  = (bext*q_mass).rename('aot')
         else:
            print('aot needs q_mass as input')

      if name == "ssa":
         bext = self.getVariable('bext', bin, wavelength = wavelength, channel=channel)           
         bsca = self.getVariable('bsca', bin, wavelength = wavelength, channel=channel)           
         ssa  = bsca/bext
         var  = ssa.interp(rh = rh).rename('ssa')

      if name == 'volume':
         rhod = self.getVariable('rhod', bin)
         gf   = self.getVariable('gf', bin)
         if rhod is not None and gf is not None:
            vol = gf**3/rhod
            var = vol.interp(rh = rh).rename('volume')
         else:
            print('vol needs rhod and growth factor in the table')

      if name == 'area':
         rhod  = self.getVariable('rhod', bin)
         gf    = self.getVariable('gf', bin)
         reff  = self.getVariable('rEff', bin)
         if (rhod is not None and gf is not None and reff is not None):
            vol  = gf**3/rhod
            area = vol/(4./3.*reff)
            var  = area.interp(rh = rh).rename('area')
         else:
            print('area needs rhod, growth factor and rEff in the table')

      if name == 'RefIndex':
         refr = self.getVariable('refreal', bin, rh = rh, wavelength = wavelength, channel = channel)
         refi = self.getVariable('refimag', bin, rh = rh, wavelength = wavelength, channel = channel)
         var = (refr, refi)

      if name == 'aot_ssa_pmom' or name == 'aot_ssa_gasym':
         if q_mass is not None:
            bext = self.getVariable('bext', bin, wavelength = wavelength, channel= channel)
            bsca = self.getVariable('bsca', bin, wavelength = wavelength, channel= channel)
            ssa  = (bsca/bext).interp(rh = rh).rename('ssa')
            aot  = (bext.interp(rh = rh) * q_mass).rename('aot')
            if 'pmom' in name:
               pmom = self.getVariable('pmom', bin, rh = rh, wavelength = wavelength, channel= channel)
               var  = (aot, ssa, pmom)
            if 'gasym' in name:
               gasym = self.getVariable('g', bin, rh = rh, wavelength = wavelength, channel= channel)
               var   = (aot, ssa, gasym)
         else:
            print( name + ' needs q_mass as input')

      if name == 'p11':
         p11 = self.getVariable('pback', bin, wavelength = wavelength, channel= channel)
         if (p11 is not None):
            var = p11.isel({"nPol": [0]}).interp(rh = rh).rename('p11')
         else:
            print('p11 needs pback in the table')

      if name == 'p22':
         p22 = self.getVariable('pback', bin, wavelength = wavelength, channel= channel)
         if (p22 is not None):
            var = p22.isel({"nPol": [4]}).interp(rh = rh).rename('p22')
         else:
            print('p22 needs pback in the table')
 
      return var
       
#--
   def getChannel(self, wavelength):
      """
        Returns channel number for a given wavelength. If it is not found, return None
      """
      tol_    = 10e-9
      channel = self.mieDS.indexes['lambda'].get_indexer([wavelength], method = 'nearest', tolerance = tol_)[0] 
      return channel if channel != -1 else None 

#--
   def getWavelength(self, channel):
      """
        Returns channel number for a given wavelength.
      """
      wavelength = this.wavelengths.values[channel]
      return wavelength

#--
   def getAOT(self, q_mass, rh, bin, wavelength=None, channel=None):
      """
        Aerol extinction optical depth (used to be called TAU)
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)
  
      bext_ = self.bext.isel({'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      AOT   = bext_*q_mass
      return AOT

#--
   def getGASYM(self, rh, bin, wavelength=None, channel=None):
      """
        Asymmetry parameter
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)

      GASYM = self.g.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      return GASYM

#--
   def getBEXT(self, rh, bin, wavelength=None, channel=None):
      """
        Mass extinction efficiency [m2 (kg dry mass)-1]
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)

      BEXT = self.bext.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      return BEXT

#--
   def getBSCA(self, rh, bin, wavelength=None, channel=None):
       """
        Mass scattering efficiency [m2 (kg dry mass)-1]
       """
       bin_     = bin - 1
       channel_ = channel
       if not channel_:
          assert wavelength, "Either wavelength or channel should be provided as input"
          channel_= self.getChannel(wavelength)

       BSCA = self.bsca.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
       return BSCA

#--
   def getSSA(self, rh, bin, wavelength=None, channel=None):
      """
       Single scattering albedo
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)

      bext_ = self.bext.isel( {'radius':[bin_], 'lambda':[channel_]})
      bsca_ = self.bsca.isel( {'radius':[bin_], 'lambda':[channel_]})
      SSA   = (bsca_/bext_).interp(rh = rh)
      return SSA

#--
   def getBBCK(self, rh, bin, wavelength=None, channel=None):
      """
        Mass backscatter efficiency [m2 (kg dry mass)-1]
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)

      BBCK  = self.bbck.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      return BBCK

#--
   def getREFF(self, rh, bin):
      """
       Effective radius (micron) 
      """
      bin_   = bin - 1
      REFF   = self.reff.isel({'radius':[bin_]}).interp(rh = rh)
      return REFF

#--
   def getPMOM(self, rh, bin, wavelength=None, channel=None):
      """
       Moments of phase function
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)

      PMOM  = self.pmom.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      return PMOM

#--
   def getGF(self, rh, bin):
      """
       Growth factor (ratio of wet to dry radius) 
      """
      bin_ = bin - 1
      GF   = self.gf.isel( {'radius':[bin_]}).interp(rh = rh)
      return GF

#--
   def getRHOP(self, rh, bin):
      """
       Wet particle density [kg m-3] 
      """
      bin_  = bin - 1
      RHOP  = self.rhop.isel( {'radius':[bin_]}).interp(rh = rh)
      return RHOP

#--
   def getRHOD(self, rh, bin):
      """
       Dry particle density [kg m-3] 
      """
      bin_  = bin - 1
      RHOD  = self.rhod.isel( {'radius':[bin_]}).interp(rh = rh)
      return RHOD
   
#--
   def getVOLUME(self, rh, bin):
      """
        Wet particle volume [m3 kg-1]
      """
      bin_   = bin - 1
      rhod   = self.rhod.isel({'radius':[bin_]})
      gf     = self.gf.isel({'radius':[bin_]})
      vol    = gf**3/rhod
      VOLUME = vol.interp(rh = rh)
      return VOLUME
   
#--
   def getAREA(self, rh, bin):
      """
        Wet particle cross section [m2 kg-1]
      """
      bin_   = bin - 1

      rhod   = self.rhod.isel({'radius':[bin_]})
      gf     = self.gf.isel({'radius':[bin_]})
      vol    = gf**3/rhod
      reff   = self.reff.isel({'radius':[bin_]})
      area   = vol /(4./3.*reff)
      AREA   = area.interp(rh = rh)
      return AREA
   
#--
   def getVECTOR(self, q_mass, rh, bin, wavelength=None, channel=None):
      """
        vector of (aot, ssa, pmom)
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "either provide wavelength or channel"
         channel_= self.getChannel(wavelength)

      bext_ = self.bext.isel( {'radius':[bin_], 'lambda':[channel_]})
      bsca_ = self.bsca.isel( {'radius':[bin_], 'lambda':[channel_]})
      ssa   = bsca_/bext_        
      AOT   = bext_.interp(rh = rh)*q_mass
      SSA   = ssa.interp(rh = rh)
      PMOM  = self.getPMOM(rh, bin, channel=channel_)
      return (AOT, SSA, PMOM)
#--
   def getScalar(self, q_mass, rh, bin, wavelength=None, channel=None):
      """
        (aot, ssa, gasym)
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "either provide wavelength or channel"
         channel_= self.getChannel(wavelength)

      bext_ = self.bext.isel( {'radius':[bin_], 'lambda':[channel_]})
      bsca_ = self.bsca.isel( {'radius':[bin_], 'lambda':[channel_]})
      AOT   = bext_.interp(rh = rh)*q_mass
      SSA   = (bsca_/bext_).interp(rh = rh)
      GASYM = self.g.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      return (AOT, SSA, GASYM)

#--
   def getRefIndex(self, rh, bin, wavelength=None, channel=None):
      """
        Mass backscatter efficiency [m2 (kg dry mass)-1]
      """
      bin_     = bin - 1
      channel_ = channel
      if not channel_:
         assert wavelength, "Either wavelength or channel should be provided as input"
         channel_= self.getChannel(wavelength)

      REFR  = self.refr.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      REFI  = self.refi.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
      return (REFR, REFI) 
   
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
   
   # Dust
   # ----
   table  = Tables[0]
   #wavelengths = [470e-9, 550e-9, 670e-9, 870e-9]
   mie    = MieTABLE(table)#,wavelengths)
   q_mass = aer['DU003'] * aer['DELP'] / 9.81
   rh     = aer['RH']

   aot    = mie.getAOT(q_mass,rh, 3, wavelength=550e-9)
   aot1   = mie.getVariable('aot', 3, rh=rh, q_mass=q_mass, wavelength=550e-9)

   gasym  = mie.getGASYM(rh, 3, wavelength=550e-9)
   gasym1 = mie.getVariable('gasym', 3, rh=rh, wavelength=550e-9)

   bext   = mie.getBEXT(rh, 3, wavelength=550e-9)
   bext1  = mie.getVariable('bext', 3, rh=rh, wavelength=550e-9)

   bsca   = mie.getBSCA(rh, 3, wavelength=550e-9)
   bsca1  = mie.getVariable('bsca', 3, rh=rh, wavelength=550e-9)

   ssa    = mie.getSSA(rh, 3, wavelength=550e-9)
   ssa1   = mie.getVariable('ssa', 3, rh=rh, wavelength=550e-9)

   bbck   = mie.getBBCK(rh, 3, wavelength=550e-9)
   bbck1  = mie.getVariable('bbck', 3, rh=rh, wavelength=550e-9)

   reff   = mie.getREFF(rh, 3)
   reff1  = mie.getVariable('rEff', 3, rh=rh)

   (refr, refi)         = mie.getRefIndex(rh, 3, wavelength=550e-9) 
   (refr1, refi1)       = mie.getVariable('RefIndex', 3, rh=rh, wavelength=550e-9)

   (aot, ssa,gasym)     = mie.getScalar(q_mass,rh, 3, wavelength=550e-9)
   (aot1, ssa1, gasym1) = mie.getVariable('aot_ssa_gasym', 3, q_mass=q_mass, rh=rh, wavelength=550e-9)

   vol  = mie.getVariable('volume', 3, rh=rh)
   p11  = mie.getVariable('p11', 3, rh=rh, wavelength=550e-9)
   p21  = mie.getVariable('p21', 3, rh=rh, wavelength=550e-9)
   #pmom  = mie.getVariable('pmom', 3, q_mass=q_mass, rh=rh, wavelength=550e-9)

   #gf    = mie.getGF(rh, 3)
   #rhop  = mie.getRHOP(rh, 3)
   #rhod  = mie.getRHOD(rh, 3)
   #pmom  = mie.getPMOM(rh, 3, wavelength=550e-9)
   #vol = mie.getVOLUME(rh,3)
   #area  = mie.getAREA(rh, 3)
   #(aot, ssa, pmom) = mie.getVECTOR(q_mass,rh,3, wavelength=550e-9)

   # OC
   # --
   #table = Tables[1]
   #wavelengths = [470e-9, 550e-9, 670e-9, 870e-9]
   #mie = MieTABLE(table)
   
   #q_mass = aer['BCPHILIC'] * aer['DELP'] / 9.81
   #rh = aer['RH']
   #aot = mie.getAOT(q_mass,rh, 2, wavelength=550e-9)
   #vol = mie.getVOLUME(rh,3)
   #(aot, ssa, pmom) = mie.getVECTOR(q_mass,rh,3, wavelength=550e-9) 
