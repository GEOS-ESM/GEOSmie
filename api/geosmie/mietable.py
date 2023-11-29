#!/usr/bin/env python3
#

"""

   Implements API to access Mie LUTs for a single specie.

"""

import numpy as np
import xarray as xr
import pandas as pd

__VERSION__ = '0.9.0'

class MieTABLE(object):

    def __init__ ( self, filename, wavelengths=None, nmom=None ):
        """
        Uses xarray to load Mie table from file name. On inout

        filename    --- Mie Table file name
        wavelengths --- desired wavelengths; if omitted all wavelengths on file
        nmon        --- number of moments for phase function; if omitted no phase
                        function needed (saves memory).

        """
        mieDataset = xr.open_dataset(filename)

        # Single Scattering properties (fix order of dimension, what is below comes from Fortran
        # ----------------------------
        self.wavelengths = mieDataset.coords['lambda']     # (c) wavelengths        [m]
        self.rh          = mieDataset.coords['rh']         # (r) RH values   [fraction]
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
        
        vars_list = mieDataset.keys()

        if 'refreal' in vars_list:
           self.refr       = mieDataset.refreal   
        if 'refimag' in vars_list:
           self.refi       = mieDataset.refimag
        if 'pmom' in vars_list:
           self.pmom       = mieDataset.pmom
        if 'bbck' in vars_list:
           self.bbck       = mieDataset.bbck
        if 'g' in vars_list:
           self.g          = mieDataset.g
        if 'bext' in vars_list:
           self.bext       = mieDataset.bext
        if 'bsca' in vars_list:
           self.bsca       = mieDataset.bsca
        if 'rhop' in vars_list:
           self.rhop       = mieDataset.rhop
        if 'rhod' in vars_list:
           self.rhod       = mieDataset.rhod
        if 'rEff' in vars_list:
           self.reff       = mieDataset.rEff
        if 'growth_factor' in vars_list:
           self.gf         = mieDataset.growth_factor

#---
    def getXXX(self, q_mass, rh, bin, wavelength=None, channel=None):
        """
        Given either channel or wavelength, compute XXX,
        preserving the shape of the input:

        q_mass ---  aerosol layer mass profile
        rh     ---  relative humidity [0,1]

        Notice that q_mass and rh must have the same shape.
        
        """

        return None
        
#---
    def getChannel(self, wavelength):
        """
        Returns channel number for a given wavelength.
        """
        t_tol   = 10e-9
        channel = 0
        for w in self.wavelengths.values:
          if (wavelength-t_tol <= w and w<=wavelength+t_tol):
            return channel
          channel +=1
        return None
#---
    def getWavelength(self, channel):
        """
        Returns channel number for a given wavelength.
        """
        wavelength = this.wavelengths.values[channel]
        return wavelength

#---
    def getAOT(self, q_mass, rh, bin, wavelength=None, channel=None):
        """
          Aerol extinction optical depth (used to be called TAU)
        """
        bin_     = bin - 1
        channel_ = channel
        if not channel_:
           assert wavelength, "Either wavelength or channel should be provided as input"
           channel_= self.getChannel(wavelength)
  
        bext_ = self.bext.isel( {'radius':[bin_], 'lambda':[channel_]}).interp(rh = rh)
        AOT   = bext_*q_mass
        return AOT

#---
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

#---
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

#---
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

#---
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

#---
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

#---
    def getREFF(self, rh, bin):
        """
         Effective radius (micron) 
        """
        bin_     = bin - 1
        REFF  = self.reff.isel({'radius':[bin_]}).interp(rh = rh)
        return REFF

#---
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

#---
    def getGF(self, rh, bin):
        """
         Growth factor (ratio of wet to dry radius) 
        """
        bin_ = bin - 1
        GF   = self.gf.isel( {'radius':[bin_]}).interp(rh = rh)
        return GF

#---
    def getRHOP(self, rh, bin):
        """
         Wet particle density [kg m-3] 
        """
        bin_  = bin - 1
        RHOP  = self.rhop.isel( {'radius':[bin_]}).interp(rh = rh)
        return RHOP

#---
    def getRHOD(self, rh, bin):
        """
         Dry particle density [kg m-3] 
        """
        bin_  = bin - 1
        RHOD  = self.rhod.isel( {'radius':[bin_]}).interp(rh = rh)
        return RHOD
    
#---
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
    
#---
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
    
#---
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
#---
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

#---
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
   gasym  = mie.getGASYM(rh, 3, wavelength=550e-9)
   bext   = mie.getBEXT(rh, 3, wavelength=550e-9)
   bsca   = mie.getBSCA(rh, 3, wavelength=550e-9)
   ssa    = mie.getSSA(rh, 3, wavelength=550e-9)
   bbck   = mie.getBBCK(rh, 3, wavelength=550e-9)
   reff   = mie.getREFF(rh, 3)
   (refr, refi)     = mie.getRefIndex(rh, 3, wavelength=550e-9) 
   (aot, ssa,gasym) = mie.getScalar(q_mass,rh, 3, wavelength=550e-9)
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
