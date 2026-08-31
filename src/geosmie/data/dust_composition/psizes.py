#!/usr/bin/env python3

# Some codes to figure out part sizes needed for GEOSmie based
# on Mian's provided parameters

import numpy as np

nbin  = 8
sigmaarr = [1.5,1.5,1.5,1.5,1.8,2.0,2.2,2.4]
rhop  = 2600.  # kg m-3
reff  = [0.14, 0.24, 0.45, 0.8, 1.4, 2.4, 4.5, 8.0]  # effective radius [um]
cpi = 4./3.*np.pi

def carmabins(nbin, rmrat, rmin, rhop=1.):

    """
     Procedure returns a carma-like distributions of radius bins
     The bins are centered in volume betwen rlow and rup
     That is:
      r^3-rlow^3 = rup^3-r^3,
     or equivalently
      r^3. = (rup^3.+rlow^3.)/2.
     rup^3 = rlow^3*rmrat, which can be solved to find r given rmrat
     and a desired rlow, e.g. rmin = 1.d-6*((1.+rmrat)/2.)^(1.d/3)
     where the desired rlow = 1.d-6 in this example.
     The meaning of r is that it is the radius of the particle with
     the average volume of the bin.

     Variables:
     Input
      nbin = number of size bins desired
      rmrat = ratio of volume (mass) between size bins
      rmin = radius (r) of smallest bins
      rhop = particle density
     Output
      rmass = mass of bin (4./3.*pi*r^3.)*rhop
      rmassup = mass of upper limit of bin (4./3.*pi*rup^3.)*rhop
      r = radius
      rup = bin upper edge radius
      dr = width of bin (rup - rlow)
      rlow = bin lower edge radius

     This code is adapted from carmabins.pro IDL code which is itself
     based on CARMA setupbins.f90
    """
    rmassmin = cpi*rhop*rmin**3.
    vrfact = ( (3./2./np.pi / (rmrat+1))**(1./3.))*(rmrat**(1./3.) - 1.)

    rmass   = np.zeros(nbin)
    rmassup = np.zeros(nbin)
    r       = np.zeros(nbin)
    rup     = np.zeros(nbin)
    dr      = np.zeros(nbin)
    rlow    = np.zeros(nbin)

    for ibin in range(nbin):
        rmass[ibin]   = rmassmin*rmrat**ibin

    rmassup = 2.*rmrat/(rmrat+1.)*rmass
    r       = (rmass/rhop/cpi)**(1./3.)
    rup     = (rmassup/rhop/cpi)**(1./3.)
    dr      = vrfact*(rmass/rhop)**(1./3.)
    rlow    = rup - dr

    return r, dr, rlow, rup, rmassup



def lognormal(r, rm, sigma, N=1.):
    """
    Fill out a set of discrete bins using the parameters of a lognormal
    distribution. This function produces dN/dr form of the function. Note
    that this version takes sigma, not S; S = ln(sigma)

    Arguments:
    r -- the radii at which to evaluate the function
    rm -- median radius of the lognormal distribution
    sigma -- width parameter of the lognormal distribution
    N -- number/scaling parameter of the distribution
    """
    dndlogr = (N/(np.log(sigma)*np.sqrt(2*np.pi)))*np.exp(-1*(np.log(2*r)-np.log(2*rm))**2/(2*np.log(sigma)**2))
    return dndlogr


def rnum(reff,sigma):
    rn = reff*np.exp(-5./2.*np.log(sigma)**2.)
    return rn

def mostmass(rn,sigma,thresh=0.999):
#   Find edge radii where 99.9% of mass is captured

#   Define large number of subbins
    nbin_ = 1000
    nb2   = int(nbin_/2)
    rlow_ = rn/100.
    rup_  = rn*100.
    rmrat = (rup_**3./rlow_**3.)**(1./nbin_)
    r, dr, rlow, rup, rmassup = carmabins(nbin_,rmrat,rlow_,rhop=rhop)

#   Now get the number size distribution
    dndlogr = lognormal(r,rn,sigma)
    masstot = np.sum(r**3. * dndlogr * dr/r)

#   Now iterate from center point
    mtest = 0.
    j = 0
    while mtest < thresh*masstot:
        j += 1
        mtest = np.sum(r[nb2-j:nb2+j]**3. * dndlogr[nb2-j:nb2+j]
                         * dr[nb2-j:nb2+j] / r[nb2-j:nb2+j])
    rmin = rlow[nb2-j]
    rmax = rup[nb2+j]

    return rmin, rmax

def printv(f,val,i,format=True):
    if format:
        if i < nbin-1:
            f.write("[%10.4e],"%val)
        else:
            f.write("[%10.4e]],\n"%val)
    else:
        if i < nbin-1:
            f.write("%d,"%val)
        else:
            f.write("%d],\n"%val)
    return


if __name__ == "__main__":
    with open("out.txt","w") as f:
        f.write('  "r0": [')
        for i in range(0,nbin):
            rn = rnum(reff[i],sigmaarr[i])  # still microns
            r0 = rn/1.e6
            printv(f,r0,i)


        f.write('  "rmax0": [')
        for i in range(0,nbin):
            rn = rnum(reff[i],sigmaarr[i])  # still microns
            rmin, rmax = mostmass(rn,sigmaarr[i])
            rmax0 = rmax/1.e6
            printv(f,rmax0,i)
   
        f.write('  "rmin0": [')
        for i in range(0,nbin):
            rn = rnum(reff[i],sigmaarr[i])  # still microns
            rmin, rmax = mostmass(rn,sigmaarr[i])
            rmin0 = rmin/1.e6
            printv(f,rmin0,i)

        f.write('  "numperdec": [')
        for i in range(0,nbin):
            printv(f,400,i,format=False)

        f.write('  "sigma": [')
        for i in range(0,nbin):
           printv(f,sigmaarr[i],i)

        f.write('  "fracs": [')
        for i in range(0,nbin):
            printv(f,1.0,i)

    f.close()
