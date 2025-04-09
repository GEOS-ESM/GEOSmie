#!/usr/bin/env python3
"""
Read an geosmie optics table
Evaluate the GSF with weights from table to recreate phase matrix
Compare the recreated phase matrix and s11, s12, s22, s33, s34, s44

P. Castellanos April 2025


"""
import os
import xarray as xr
import numpy as np
from math import sqrt, cos, factorial
import matplotlib.pyplot as plt
from scipy.special import legendre


def recur_d_mn(m,n,s,theta):
    smin = max(abs(m),abs(n))
    if s < smin:
        return 0
    elif s == smin:
        if n >= m:
            xi = 1
        else:
            xi = (-1)**(m-n)
        sqrt_term = factorial(2*smin)/(factorial(abs(m-n))*factorial(abs(m+n)))
        third_term  = (1-cos(theta))**(0.5*abs(m-n))
        fourth_term = (1+cos(theta))**(0.5*abs(m+n))

        return xi*2**(-1*smin)*sqrt(sqrt_term)*third_term*fourth_term

    else:
        first_term = s*sqrt((s+1)**2 - m**2)*sqrt((s+1)**2 - n**2)
        first_term = 1/first_term
        second_term = (2*s+1)*(s*(s+1)*cos(theta) - m*n)*d_mn(m,n,s-1,theta)
        third_term = (s+1)*sqrt(s**2 - m**2)*sqrt(s**2 - n**2)*d_mn(m,n,s-2,theta)

        return first_term*(second_term - third_term)


def d_mn(m,n,s,theta,dm1=None,dm2=None):
    smin = max(abs(m),abs(n))
    if s < smin:
        return 0, 0
    elif (s == smin):
        if n >= m:
            xi = 1
        else:
            xi = (-1)**(m-n)
        sqrt_term = factorial(2*smin)/(factorial(abs(m-n))*factorial(abs(m+n)))
        third_term  = (1-cos(theta))**(0.5*abs(m-n))
        fourth_term = (1+cos(theta))**(0.5*abs(m+n))

        d = xi*(2**(-1*smin))*sqrt(sqrt_term)*third_term*fourth_term

        if (s == smin):
            return d, 0
        else:
            return d, dm1

    else:
        s = s-1
        first_term = s*sqrt((s+1)**2 - m**2)*sqrt((s+1)**2 - n**2)
        first_term = 1./first_term
        second_term = (2*s+1)*(s*(s+1)*cos(theta) - m*n)*dm1
        third_term = (s+1)*sqrt(s**2 - m**2)*sqrt(s**2 - n**2)*dm2

        d = first_term*(second_term - third_term)

        return d, dm1

class OPTICS(object):
    def __init__(self,inFile):
        """
        read in optics file and regenerate the phase matrix elements
        inFile : GEOSmie generated optics file
        """

        self.optics = xr.open_dataset(inFile)
        self.nMom = self.optics.sizes['m']
        self.theta = np.radians(self.optics['ang'])
        self.angle = self.optics['ang'].values
        self.getPmatrix()
        self.fname = os.basename(inFile)

    def getPmatrix(self):
        """
        regenerate the phase matrix from the GSF coefficients
        """

        # calculate P11
        ipol = np.argmin(np.abs(self.optics.p.values - 11))
        self.p11 = self.calcP11(ipol)

        # calculate P12
        ipol = np.argmin(np.abs(self.optics.p.values - 12))
        self.p12 = self.calcP12(ipol)

        # calculate P44
        # this is also legendre expansion by definition,
        # so use same function as P11
        ipol = np.argmin(np.abs(self.optics.p.values - 44))
        self.p44 = self.calcP11(ipol)

    def calcP11(self,ipol):
        """
        calculate P11 
        also works for P44
        """
        
        mu    = np.cos(self.theta)
        p11   = xr.DataArray(dims=self.optics['s11'].dims,coords=self.optics['s11'].coords)
        p11[:] = 0.0
        for s in range(self.nMom):
            P = legendre(s)            
            p11  += self.optics['pmom'].isel(p=ipol,m=s)*P(mu)

        return p11

    def calcP12(self,ipol):
        """
        calculate P12
        """
        m = 0
        n = 2
        p12 = xr.DataArray(dims=self.optics['s12'].dims,coords=self.optics['s12'].coords)
        p12[:] = 0.0
        for i,t in enumerate(self.theta):
            d0, dneg1 = d_mn(m,n,0,t)
            pfunc = d0/(1j**(n-m)).real
            p12.isel(ang=i).values = self.optics['pmom'].isel(m=0,p=ipol)*pfunc
            d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
            pfunc = d1/(1j**(n-m)).real
            p12.isel(ang=i).values += self.optics['pmom'].isel(m=1,p=ipol)*pfunc

            dm2 = d0
            dm1 = d1
            for s in range(2,self.nMom):
                dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                pfunc = dfunc/(1j**(n-m)).real

                p12.isel(ang=i).values += self.optics['pmom'].isel(m=s,p=ipol)*pfunc

                dm1 = dfunc

        return p12
        
    def plotPmatrix(self,irh=None,iwav=None,ibin=None):
        """
        plot the phase matrices
        can optionally provide a single rh, wavelength, or bin number to plot
        otherwise will plot all
        """

        if irh is None:
            irh = np.arange(self.optics.rh.size)
        if iwav is None:
            iwav = np.arange(self.optics.wavelength.size)
        if ibin is None:
            ibin = np.arange(self.optics.bin.size)

        for rh in irh:
            for wav in iwav:
                for bin in ibin:
                    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))
                    # P11
                    ax = axes[0,0]
                    p11 = self.p11.isel(rh=rh,wavelength=wav,bin=bin)    
                    ax.plot(self.angle, p11,label='GSF')
                    # normalize S11
                    s11 = self.optics['s11'].isel(rh=rh,wavelength=wav,bin=bin)
                    s11n = 2.*s11 / np.trapz(s11 * np.sin(self.theta),self.theta)
                    ax.plot(self.angle,s11n,label='S11')
                    ax.legend()
                    ax.set_title('P11=P1')           

                    # P12
                    ax = axes[0,1]
                    ax.plot(self.angle, self.p12.isel(rh=rh,wavelength=wav,bin=bin)/p11,label='GSF')
                    ax.plot(self.angle,self.optics['s12'].isel(rh=rh,wavelength=wav,bin=bin)/s11,label='S12')
                    ax.legend()
                    ax.set_title('P12=P2')

                    plt.show()

 

if __name__ == '__main__':

    inDir =  '/discover/nobackup/pcastell/workspace/aero_work/aist/sbg/aop_testing/ExtDataColarco'
    inFile = 'optics_SU.v2.0.0.GSFun-129.nc4'


    optics = OPTICS(inDir+'/'+inFile)
    optics.plotPmatrix(irh=[0],ibin=[0],iwav=[0])
