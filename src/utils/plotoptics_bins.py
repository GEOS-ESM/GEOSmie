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
import matplotlib
#matplotlib.use('agg')
from matplotlib.ticker import ScalarFormatter

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def singleplot(ax,x,y,ylabel,title,i40,i80,i95,waveum,legend=False):
        def diamtosize(d):
            return np.pi*d/waveum
        def sizetodiam(x):
            return waveum/np.pi/x
        ax.set_xlabel(r'effective diameter [$\mu$m]')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if(legend):
            ax.plot(x[:,0],   y[:,0],label="RH=0%")
            ax.plot(x[:,i40], y[:,i40],label="RH=40%")
            ax.plot(x[:,i80], y[:,i80],label="RH=80%")
            ax.plot(x[:,i95], y[:,i95],label="RH=95%")
            ax.legend()
        else:
            ax.plot(x[:,0],   y[:,0])
            ax.plot(x[:,i40], y[:,i40])
            ax.plot(x[:,i80], y[:,i80])
            ax.plot(x[:,i95], y[:,i95])
        ax.set_xscale('log')
        ax2 = ax.secondary_xaxis("top",functions=(diamtosize,sizetodiam))
        ax2.set_xlabel('Size Parameter')
        ax2.set_xscale('log')
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax2.xaxis.set_major_formatter(ScalarFormatter())

class OPTICS(object):
    def __init__(self,inFile,irh=None,ibin=None,iwav=None):
        """
        read in optics file and regenerate the phase matrix elements
        inFile : GEOSmie generated optics file


        irh, ibin, iwav = limit calculations to these indices, if None do all
        """

        self.fname = os.path.basename(inFile)
        
        # read optics file
        self.optics = xr.load_dataset(inFile,engine="netcdf4")
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

    def plot(self, wavelength=550, ibin=None,fig=None,axes=None,plotclose=True):
        """
        plot the single scattering properties
        assume all RH on file
        can optionally provide a single bin number to plot
        otherwise will plot all
        """

        print(self.optics['bext'].shape)
        iwav = find_nearest(self.wavelengths,wavelength)
        i450 = find_nearest(self.wavelengths,450)
        i900 = find_nearest(self.wavelengths,900)
        wave = self.wavelengths[iwav]
        w450 = self.wavelengths[i450]
        w900 = self.wavelengths[i900]
        
        i40 = find_nearest(self.rh,0.4)
        i80 = find_nearest(self.rh,0.8)
        i95 = find_nearest(self.rh,0.95)
        
        waveum = wave/1000.
        def diamtosize(d):
            return np.pi*d/waveum
        def sizetodiam(x):
            return waveum/np.pi/x

        if ibin is None:
            ibin = np.arange(self.optics['bext'].sizes['bin'])

        if(fig==None):
            fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(16, 12))

        bext = self.optics['bext'].isel(wavelength=iwav)
        bsca = self.optics['bsca'].isel(wavelength=iwav)
        bbck = self.optics['bbck'].isel(wavelength=iwav)
        b450 = self.optics['bext'].isel(wavelength=i450)
        b900 = self.optics['bext'].isel(wavelength=i900)
        g    = self.optics['g'].isel(wavelength=iwav)
        reff = self.optics["rEff"]*1e6
        scatang = self.optics["ang"]
        p11  = self.optics["p11"].isel(wavelength=iwav)
        p12  = self.optics["p12"].isel(wavelength=iwav)
        ang  = -np.log(b450/b900)/np.log(w450/w900)
        sr   = self.optics["lidar_ratio"].isel(wavelength=iwav)

        # BEXT as function of size
        ax = axes[0,0]
        x = 2.*reff
        singleplot(ax,x,bext,r'$\beta_{ext} [m^2 g^{-1}]$',
                  'Extinction Efficiency [%3d nm]'%(wave),i40,i80,i95,waveum,legend=True)
        
        # SSA as function of size
        ax = axes[0,1]
        singleplot(ax,x,bsca/bext,"SSA",
                   'Single Scattering Albedo [%3d nm]'%(wave),i40,i80,i95,waveum)

        # Lidar ratio
        ax = axes[0,2]
        singleplot(ax,x,sr,r'$s_{r} [sr^{-1}]$',
                   'Lidar Ratio [%3d nm]'%(wave),i40,i80,i95,waveum)

        # Asymmetry factor as function of size
        ax = axes[1,0]
        singleplot(ax,x,g,"g",
                   'Asymmetry Parameter [%3d nm]'%(wave),i40,i80,i95,waveum)

        # P11
        ax = axes[1,1]
        titlestr = 'Phase Function [%3d nm]'%(wave)
        ax.set_xlabel('scattering angle')
        ax.set_ylabel('P11')
        ax.set_title(titlestr)
        for jbin in ibin:
            ax.plot(scatang,p11[jbin,0,:],label=f"bin {jbin}")
        ax.set_yscale('log')
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.legend()


        # DOLP
        ax = axes[1,2]
        titlestr = 'Degree of Linear Polarization [%3d nm]'%(wave)
        ax.set_xlabel('scattering angle')
        ax.set_ylabel('DOLP')
        ax.set_title(titlestr)
        for jbin in ibin:
            ax.plot(scatang,-p12[jbin,0,:]/p11[jbin,0,:],label=f"bin {jbin}")
        ax.legend()


       
        plt.suptitle(self.fname+' [%3d nm]'%(wave))
        plt.tight_layout(pad=4)
        if(plotclose):
            plt.savefig('plots/' + self.fname + '_wav%04dnm.png'%(wave))
        return fig, axes


if __name__ == '__main__':
    parser = OptionParser(usage="Usage: %prog",
                          version='0.0.1' )
    parser.add_option("--name", dest="name", default="",
                      help="Particle file to plot (default=%s)"\
                      %(""))
    (options, args) = parser.parse_args()
    fname = options.name
  
    optics = OPTICS(fname)
    fig, axes = optics.plot(wavelength=10000)
    optics.plot()

