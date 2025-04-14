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
    def __init__(self,inFile,irh=None,ibin=None,iwav=None):
        """
        read in optics file and regenerate the phase matrix elements
        inFile : GEOSmie generated optics file


        irh, ibin, iwav = limit calculations to these indices, if None do all
        """

        self.fname = os.path.basename(inFile)
        
        # read optics file
        self.optics = xr.load_dataset(inFile)
        self.nMom = self.optics.sizes['m']
        self.theta = np.radians(self.optics['ang'])
        self.angle = self.optics['ang'].values

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
                       'ang': self.optics['ang'],
                      } 

        # calculate phase matrix from expansion coefficients
        self.getPmatrix()

        # pull out the indices of S-elements selected for plotting later
        self.s11 = self.optics['s11'].isel(rh=self.irh,bin=self.ibin,wavelength=self.iwav)
        self.s12 = self.optics['s12'].isel(rh=self.irh,bin=self.ibin,wavelength=self.iwav)
        self.s22 = self.optics['s22'].isel(rh=self.irh,bin=self.ibin,wavelength=self.iwav)
        self.s33 = self.optics['s33'].isel(rh=self.irh,bin=self.ibin,wavelength=self.iwav)
        self.s34 = self.optics['s34'].isel(rh=self.irh,bin=self.ibin,wavelength=self.iwav)
        self.s44 = self.optics['s44'].isel(rh=self.irh,bin=self.ibin,wavelength=self.iwav)

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

        # calculate P22 & P33
        ipol_22 = np.argmin(np.abs(self.optics.p.values - 22))
        ipol_33 = np.argmin(np.abs(self.optics.p.values - 33))
        self.p22, self.p33 = self.calcP22_P33(ipol_22,ipol_33)

        # calculate P34
        ipol = np.argmin(np.abs(self.optics.p.values - 34))
        self.p34 = self.calcP34(ipol)

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
        p11   = xr.DataArray(dims=self.dims,coords=self.coords)
        p11[:] = 0.0
        for s in range(self.nMom):
            P = legendre(s)            
            p11  += self.optics['pmom'].isel(p=ipol,m=s,rh=self.irh,wavelength=self.iwav,bin=self.ibin)*P(mu)

        return p11

    def calcP12(self,ipol):
        """
        calculate P12
        """
        m = 0
        n = 2
        p12 = xr.DataArray(dims=self.dims,coords=self.coords)
        p12[:] = 0.0
        for i,t in enumerate(self.theta):
            d0, dneg1 = d_mn(m,n,0,t)
            pfunc = d0/(1j**(n-m)).real
            p12[dict(ang=i)] = self.optics['pmom'].isel(m=0,p=ipol,rh=self.irh,wavelength=self.iwav,bin=self.ibin)*pfunc
            d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
            pfunc = d1/(1j**(n-m)).real
            p12[dict(ang=i)] += self.optics['pmom'].isel(m=1,p=ipol,rh=self.irh,wavelength=self.iwav,bin=self.ibin)*pfunc

            dm2 = d0
            dm1 = d1
            for s in range(2,self.nMom):
                dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                pfunc = dfunc/(1j**(n-m)).real

                p12[dict(ang=i)] += self.optics['pmom'].isel(m=s,p=ipol,rh=self.irh,wavelength=self.iwav,bin=self.ibin)*pfunc

                dm1 = dfunc

        return p12

    def calcP22_P33(self,ipol_22,ipol_33):
        """
        calculate P22 & P33
        """

        # first get a2 + a3
        m = 2
        n = 2
        a2p3 = xr.DataArray(dims=self.dims,coords=self.coords)
        a2p3[:] = 0.0
        gsf_coef22 = self.optics['pmom'].isel(p=ipol_22,rh=self.irh,wavelength=self.iwav,bin=self.ibin)
        gsf_coef33 = self.optics['pmom'].isel(p=ipol_33,rh=self.irh,wavelength=self.iwav,bin=self.ibin)
        for i,t in enumerate(self.theta):
            d0, dneg1 = d_mn(m,n,0,t)
            pfunc = d0/(1j**(n-m)).real
            a2p3[dict(ang=i)] = (gsf_coef22.isel(m=0)+gsf_coef33.isel(m=0))*pfunc
            d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
            pfunc = d1/(1j**(n-m)).real
            a2p3[dict(ang=1)] += (gsf_coef22.isel(m=1)+gsf_coef33.isel(m=1))*pfunc

            dm2 = d0
            dm1 = d1
            for s in range(2,self.nMom):
                dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                pfunc = dfunc/(1j**(n-m)).real

                a2p3[dict(ang=i)] +=  (gsf_coef22.isel(m=s)+gsf_coef33.isel(m=s))*pfunc

                dm1 = dfunc

        # next get a2 - a3
        m = 2
        n = -2
        a2m3 = xr.DataArray(dims=self.dims,coords=self.coords)
        a2m3[:] = 0.0
        for i,t in enumerate(self.theta):
            d0, dneg1 = d_mn(m,n,0,t)
            pfunc = d0/(1j**(n-m)).real
            a2m3[dict(ang=i)] = (gsf_coef22.isel(m=0)-gsf_coef33.isel(m=0))*pfunc
            d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
            pfunc = d1/(1j**(n-m)).real
            a2m3[dict(ang=i)] += (gsf_coef22.isel(m=1)-gsf_coef33.isel(m=1))*pfunc

            dm2 = d0
            dm1 = d1
            for s in range(2,self.nMom):
                dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                pfunc = dfunc/(1j**(n-m)).real

                a2m3[dict(ang=i)] += (gsf_coef22.isel(m=s)-gsf_coef33.isel(m=s))*pfunc

                dm1 = dfunc

        # now combine a2 + a3 = a2p3 with a2 - a3 = a2m3
        p22 = 0.5*(a2p3 + a2m3)
        p33 = a2p3 - p22

        return p22,p33

    def calcP34(self,ipol):
        """
        calculate P34
        """

        m = 0
        n = 2

        p34 = xr.DataArray(dims=self.dims,coords=self.coords)
        p34[:] = 0.0
        gsf_coef = self.optics['pmom'].isel(p=ipol,rh=self.irh,wavelength=self.iwav,bin=self.ibin)
        for i,t in enumerate(self.theta):
            d0, dneg1 = d_mn(m,n,0,t)
            pfunc = d0/(1j**(n-m)).real
            p34[dict(ang=i)] = gsf_coef.isel(m=0)*pfunc
            d1, d0 = d_mn(m,n,1,t,dm1=d0,dm2=dneg1)
            pfunc = d1/(1j**(n-m)).real
            p34[dict(ang=i)] +=  gsf_coef.isel(m=1)*pfunc

            dm2 = d0
            dm1 = d1
            for s in range(2,self.nMom):
                dfunc, dm2 = d_mn(m,n,s,t,dm1=dm1,dm2=dm2)
                pfunc = dfunc/(1j**(n-m)).real

                p34[dict(ang=i)] += gsf_coef.isel(m=s)*pfunc

                dm1 = dfunc

        return p34
        
    def plotPmatrix(self,irh=None,iwav=None,ibin=None):
        """
        plot the phase matrices
        can optionally provide a single rh, wavelength, or bin number to plot
        otherwise will plot all
        """

        if irh is None:
            irh = np.arange(self.p11.sizes['rh'])
        if iwav is None:
            iwav = np.arange(self.p11.sizes['wavelength'])
        if ibin is None:
            ibin = np.arange(self.p11.sizes['bin'])

        for rh in irh:
            for wav in iwav:
                for bin in ibin:
                    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(8, 6))
                    # P11
                    ax = axes[0,0]
                    p11 = self.p11.isel(rh=rh,wavelength=wav,bin=bin)    
                    ax.plot(self.angle, p11,label='GSF')
                    # normalize S11
                    s11 = self.s11.isel(rh=rh,wavelength=wav,bin=bin)
                    s11n = 2.*s11 / np.trapz(s11 * np.sin(self.theta),self.theta)
                    ax.plot(self.angle,s11n,label='S11')
                    ax.legend()
                    ax.set_title('P11=P1')           

                    # P12
                    ax = axes[0,1]
                    ax.plot(self.angle, self.p12.isel(rh=rh,wavelength=wav,bin=bin)/p11,label='GSF')
                    ax.plot(self.angle,self.s12.isel(rh=rh,wavelength=wav,bin=bin)/s11,label='S12')
                    ax.legend()
                    ax.set_title('P12=P2')

                    # P22
                    ax = axes[1,1]
                    ax.plot(self.angle, self.p22.isel(rh=rh,wavelength=wav,bin=bin)/p11,label='GSF')
                    ax.plot(self.angle,self.s22.isel(rh=rh,wavelength=wav,bin=bin)/s11,label='S22')
                    ax.legend()
                    ax.set_yscale('log')
                    ax.set_title('P22=P5')

                    # P33
                    ax = axes[0,2]
                    ax.plot(self.angle, self.p33.isel(rh=rh,wavelength=wav,bin=bin)/p11,label='GSF')
                    ax.plot(self.angle,self.s33.isel(rh=rh,wavelength=wav,bin=bin)/s11,label='S33')
                    ax.legend()
                    ax.set_title('P33=P3')

                    # P34
                    ax = axes[1,0]
                    ax.plot(self.angle, self.p34.isel(rh=rh,wavelength=wav,bin=bin)/p11,label='GSF')
                    ax.plot(self.angle,self.s34.isel(rh=rh,wavelength=wav,bin=bin)/s11,label='S34')
                    ax.legend()
                    ax.set_title('P34=P4')

                    # P44
                    ax = axes[1,2]
                    ax.plot(self.angle, self.p44.isel(rh=rh,wavelength=wav,bin=bin)/p11,label='GSF')
                    ax.plot(self.angle,self.s44.isel(rh=rh,wavelength=wav,bin=bin)/s11,label='S44')
                    ax.legend()
                    ax.set_title('P44=P6')

                    # Adjust layout and display
                    plt.suptitle(self.fname + ' ibin{:02d} irh={:02d} iwav={:02d}'.format(bin,rh,wav))
                    plt.subplots_adjust(wspace=0.2)
                    plt.savefig(self.fname + '_ibin{:02d}_irh{:02d}_iwav{:02d}.png'.format(bin,rh,wav))
                    #plt.show()
                    plt.close()
if __name__ == '__main__':

    inDir =  '/discover/nobackup/pcastell/workspace/aero_work/aist/sbg/aop_testing/ExtDataColarco'
#    inFile = 'optics_SU.v2.0.0.GSFun-129.nc4'
    inFile = 'optics_BC.v2.0.0.GSFun-129.nc4'


#    optics = OPTICS(inDir+'/'+inFile,irh=[0],ibin=[0],iwav=[0])
#    optics.plotPmatrix()

    optics = OPTICS(inDir+'/'+inFile,irh=[35],iwav=[0])
    optics.plotPmatrix()
