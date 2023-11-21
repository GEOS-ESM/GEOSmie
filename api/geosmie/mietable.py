"""

   Implements API to access Mie LUTs for a single specie.

"""

import numpy as np
import xarray as xr

__VERSION__ = 0.9.0

NRH_BINS = 991 

class MieTABLE(object):

    def __init__ ( self, filename, wavelengths=None, nmom=None )
        """
        Uses xarray to load Mie table from file name. On inout

        filename    --- Mie Table file name
        wavelengths --- desired wavelengths; if omitted all wavelengths on file
        nmon        --- number of moments for phase function; if omitted no phase
                        function needed (saves memory).

        """

        # Single Scattering properties (fix order of dimension, what is below comes from Fortran
        # ----------------------------
        self.wavelengths = None  # (c) wavelengths [m]
        self.rh   = None         # (r) RH values   [fraction]
        self.reff = None         # (r,b) effective radius [m]
        self.bext = None         # (r,c,b) bext values [m2 kg-1]
        self.bsca = None         # (r,c,b) bsca values [m2 kg-1]
        self.bbck = None         # (r,c,b) bbck values [m2 kg-1]
        self.g    = None         # (r,c,b) asymmetry parameter
        self.p11  = None         # (r,c,b) Backscatter phase function, index 1 
        self.p22  = None         # (r,c,b) Backscatter phase function, index 5
        self.pmom = None         # (c,r,b,m,p) moments of phase function              CHECK FORTRAN!!!!
        self.pback = None        # (c,r,b,p) moments of backscatter phase function
        self.gf   = None         # (r,b) hygroscopic growth factor
        self.rhop = None         # (r,b) wet particle density [kg m-3]
        self.rhod = None         # (r,b) dry particle density [kg m-3]
        self.vol  = None         # (r,b) wet particle volume [m3 kg-1]
        self.area = None         # (r,b) wet particle cross section [m2 kg-1]
        self.refr = None         # (r,c,b) real part of refractive index
        self.refi = None         # (r,c,b) imaginary part of refractive index

        # May or may not need this if using xarray for RH interpolation
        # -------------------------------------------------------------
        self.rhi = None      # pointer to rh LUT
        self.rha = None       # slope on rh LUT

#---
        def getChannel(self, wavelength):
            """
            Returns channel number for a given wavelength.
            """
            channel = 0
            return channel
#---
        def getWavelength(self, channel):
            """
            Returns channel number for a given wavelength.
            """
            wavelength = 0
            return wavelength

#---

        def getXXX ( self, q_mass, rh, bin, wavelength=None, channel=None):
            """
            Given either channel or wavelength, compute XXX,
            preserving the shape of the input:

            q_mass ---  aerosol layer mass profile
            rh     ---  relative humidity [0,1]

            Notice that q_mass and rh must have the same shape.
            
            """

            return

        """
        Implement getXXX for these
                          getAOT     (AOT used to be called TAU)
                          getSSA        
                          getGASYM       
                          getScalar --> AOT, SSA, GASYM
                          getBEXT
                          getBSCA
                          getBBCK
                          getREFF
                          getPMOM
                          getVector  --> AOT, SSA, PMOM
                          getGF
                          getRHOP
                          getRHOD
                          getVOLUME
                          getAREA
                          getRefIndex  --> REFR, REFI
        """
        
    #______________________________________________________________

    if __name__ == "__main__":

        # Sample Mie Tables
        # -----------------
        dirn = '/discover/nobackup/projects/gmao/share/dasilva/fvInput/ExtData/chemistry/AerosolOptics/v0.0.0/x/'
        Tables = [dirn + 'optics_DU.v7.nc', dirn + 'optics_OC.v2_3.nc']
        

        # Aerosol state (all species)
        # ---------------------------

        aer_Nv = '/css/gmao/geos-it/products/Y2023/M02/D05/GEOS.it.asm.aer_inst_3hr_glo_C180x180x6_v72.GEOS5294.2023-02-05T1200.V01.nc4'

        aer = xr.open_dataset(aer_Nv).variables
        
        # Dust
        # ----
        table = dirn + 'optics_DU.v7.nc'
        wavelengths = [470e-9, 550e-9, 670e-9, 870e-9]
        mie = MieTABLE(table,wavelengths)
        
        q_mass = aer['DU003'] * aer['DELP'] / 9.81
        rh = aer['RH']
        aot = mie.getAOT(q_mass,rh,3, 550e-9)
        vol = mie.GetVOLUME(q_mass,rh,3)
        (aot, ssa, pmom) = mie.getVECTOR(q_mass,rh,3, 550e-9)

        # OC
        # --
        table = dirn + 'optics_DU.v7.nc'
        wavelengths = [470e-9, 550e-9, 670e-9, 870e-9]
        mie = MieTABLES(table,wavelengths)
        
        q_mass = aer['BCPHILIC'] * aer['DELP'] / 9.81
        rh = aer['RH']
        aot = mie.getAOT(q_mass,rh,2, 550e-9)
        vol = mie.GetVOLUME(q_mass,rh,3)
        (aot, ssa, pmom) = mie.getVECTOR(q_mass,rh,3, 550e-9)
  
        
