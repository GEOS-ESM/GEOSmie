#!/usr/bin/env python3
"""
Read a geosmie optics table
Plot basic single scattering properties

Peter Colarco September 2025

"""
import os
import xarray as xr
import numpy as np
from math import sqrt, cos, factorial
import matplotlib.pyplot as plt
from scipy.special import eval_legendre
from optparse import OptionParser

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

class OPTICS(object):
    def __init__(self,inFile,irh=None,ibin=None,iwav=None):
        """
        read in optics file and regenerate the phase matrix elements
        inFile : GEOSmie generated optics file


        irh, ibin, iwav = limit calculations to these indices, if None do all
        """

        self.fname = os.path.basename(inFile)
        
        # read optics file
        self.optics = xr.load_dataset(inFile)
        self.rh = self.optics['rh'].values
        self.wavelengths = self.optics['wavelength'].values*1e9 # nm

        if irh is None:
            self.irh = np.arange(self.optics.sizes['rh'])
        else:
            self.irh = irh

        if ibin is None:
            self.ibin = np.arange(self.optics.sizes['bin'])
        else:
            self.ibin = ibin

        if iwav is None:
            self.iwav = np.arange(self.optics.sizes['wavelength'])
        else:
            self.iwav = iwav


        # coordinate of output
        self.dims = ('bin', 'wavelength', 'rh', 'ang')
        self.coords = {'rh': self.optics['rh'][self.irh],
                       'wavelength': self.optics['wavelength'][self.iwav],
                       'bin': self.optics['bin'][self.ibin],
                      } 

    def plot(self, wavelength=550, ibin=None):
        """
        plot the single scattering properties
        assume all RH on file
        can optionally provide a single bin number to plot
        otherwise will plot all
        """

        iwav = find_nearest(self.wavelengths,wavelength)
        i450 = find_nearest(self.wavelengths,450)
        i900 = find_nearest(self.wavelengths,900)
        wave = self.wavelengths[iwav]
        w450 = self.wavelengths[i450]
        w900 = self.wavelengths[i900]
        
        if ibin is None:
            ibin = np.arange(self.optics['bext'].sizes['bin'])

        for bin in ibin:
            fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(12, 16))

            bext = self.optics['bext'].isel(wavelength=iwav,bin=bin)
            bsca = self.optics['bsca'].isel(wavelength=iwav,bin=bin)
            bbck = self.optics['bbck'].isel(wavelength=iwav,bin=bin)
            b450 = self.optics['bext'].isel(wavelength=i450,bin=bin)
            b900 = self.optics['bext'].isel(wavelength=i900,bin=bin)
            ang  = -np.log(b450/b900)/np.log(w450/w900)
            
            # BEXT as function of RH
            ax = axes[0,0]
            ax.set_xlabel('Relative Humidity')
            ax.set_ylabel(r'$\beta_{ext} [m^2 g^{-1}]$')
            ax.set_title('Extinction Efficiency [%3d nm]'%(wave))
            ax.plot(self.rh, bext/1000.)

            # SSA as function of RH
            ax = axes[0,1]
            ax.set_xlabel('Relative Humidity')
            ax.set_ylabel('SSA')
            ax.set_title('Single Scattering Albedo [%3d nm]'%(wave))
            ax.plot(self.rh, bsca/bext)

            # Backscatter efficiency as function of RH
            ax = axes[1,0]
            ax.set_xlabel('Relative Humidity')
            ax.set_ylabel(r'$\beta_{bck} [units]$')
            ax.set_title('Backscattering Efficiency [%3d nm]'%(wave))
            ax.plot(self.rh, bbck/1000.)

            # Angstrom coefficient as function of RH
            ax = axes[1,1]
            titlestr = 'Angstrom Exponent (%3d-%3d nm pair)'%(w450,w900)
            ax.set_xlabel('Relative Humidity')
            ax.set_ylabel(r'$\alpha_{%3d-%3d}$'%(w450,w900))
            ax.set_title(titlestr)
            ax.plot(self.rh, ang)

            # BEXT as function of wavelength
            ax = axes[2,0]
            ax.set_xlim(250,5000)
            ax.set_xscale('log')
            ax.set_xlabel('Wavelength [nm]')
            ax.set_ylabel(r'$\beta_{ext} [m^2 g^{-1}]$')
            ax.set_title('Extinction Efficiency')
            i40 = find_nearest(self.rh,0.4)
            i80 = find_nearest(self.rh,0.8)
            i95 = find_nearest(self.rh,0.95)
            bext = self.optics['bext'].isel(rh=0,bin=bin)
            bsca = self.optics['bsca'].isel(rh=0,bin=bin)
            ssa00 = bsca/bext
            ax.plot(self.wavelengths, bext/1000., label='RH=0%')
            bext = self.optics['bext'].isel(rh=i40,bin=bin)
            bsca = self.optics['bsca'].isel(rh=i40,bin=bin)
            ssa40 = bsca/bext
            ax.plot(self.wavelengths, bext/1000., label='RH=40%')
            bext = self.optics['bext'].isel(rh=i80,bin=bin)
            bsca = self.optics['bsca'].isel(rh=i80,bin=bin)
            ssa80 = bsca/bext
            ax.plot(self.wavelengths, bext/1000., label='RH=80%')
            bext = self.optics['bext'].isel(rh=i95,bin=bin)
            bsca = self.optics['bsca'].isel(rh=i95,bin=bin)
            ssa95 = bsca/bext
            ax.plot(self.wavelengths, bext/1000., label='RH=95%')
            ax.legend()

            # SSA as function of wavelength
            ax = axes[2,1]
            ax.set_xlim(250,5000)
            ax.set_xscale('log')
            ax.set_title('Single Scattering Albedo')
            ax.set_xlabel('Wavelength [nm]')
            ax.set_ylabel('SSA')
            ax.plot(self.wavelengths, ssa00, label='RH=0%')
            ax.plot(self.wavelengths, ssa40, label='RH=40%')
            ax.plot(self.wavelengths, ssa80, label='RH=80%')
            ax.plot(self.wavelengths, ssa95, label='RH=95%')
            ax.legend()

            
            plt.suptitle(self.fname+' [%3d nm]'%(wave))
            plt.tight_layout(pad=4)
#            plt.show()
            plt.savefig('plots/'+self.fname + '_ibin%02d_wav%04dnm.png'%(bin,wave))

if __name__ == '__main__':
    parser = OptionParser(usage="Usage: %prog",
                          version='0.0.1' )
    parser.add_option("--name", dest="name", default="",
                      help="Particle file to plot (default=%s)"\
                      %(""))
    (options, args) = parser.parse_args()
    fname = options.name
  
    optics = OPTICS(fname)
    optics.plot()

