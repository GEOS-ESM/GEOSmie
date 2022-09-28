"""
Copyright (C) 2012-2013 Jussi Leinonen
Copyright (C) 2019-2021 Osku Kemppinen

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""

from numpy import array, arange, dot, zeros, vstack, sqrt, sin, cos
import numba
import numpy as np
import time

@numba.jit(nopython=True)
# This is still very slow and needs to be accelerated
def mie_props_raw(coeffs_raw, y):
  an, bn, nmax = coeffs_raw
  """The scattering properties.
  """
  anp = an.real
  anpp = an.imag
  bnp = bn.real
  bnpp = bn.imag

  n1 = nmax-1
  n = arange(1.,nmax+1.)
  cn = 2*n+1
  c1n = n*(n+2)/(n+1)
  c2n = cn/(n*(n+1))
  y2 = y**2

  dn = cn*(anp+bnp)
  q = dn.sum()
  qext = 2*q/y2

  en = cn*(anp**2+anpp**2+bnp**2+bnpp**2)
  q = en.sum()
  qsca = 2*q/y2
  qabs = qext-qsca

  fn = (an-bn)*cn
  gn=(-1)**n # this line takes third of the whole function exec time
  q = (fn*gn).sum()
  qb1 = q.real ** 2 + q.imag ** 2
  qb = qb1.real/y2

  g1 = zeros((4,nmax))
  g1[:,:n1] = vstack((anp[1:nmax], anpp[1:nmax], bnp[1:nmax], bnpp[1:nmax]))

  asy1 = c1n*(anp*g1[0,:]+anpp*g1[1,:]+bnp*g1[2,:]+bnpp*g1[3,:])
  asy2 = c2n*(anp*bnp+anpp*bnpp)

  asy = 4/y2*(asy1+asy2).sum()/qsca
  qratio = qb/qsca

  return (qext, qsca, qabs, qb, asy, qratio) 

def mie_props(coeffs,y):
  qext, qsca, qabs, qb, asy, qratio = mie_props_raw((coeffs.an, coeffs.bn, coeffs.nmax), y)
  return {"qext":qext, "qsca":qsca, "qabs":qabs, "qb":qb, "asy":asy, 
      "qratio":qratio}

# issues with signature, so we just use automatic signature
# we use homemade dot so we can call it from within numba-wrapped function
# fastmath improves runspeed by about 40%, which is quite meaningful (not sure about the accuracy though...)
@numba.jit(nopython=True,fastmath=False)
def numba_dot(a, b):
  ret = 0. + 0j
  for i in range(len(a)):
    ret += a[i] * b[i]
  return ret

@numba.jit(nopython=True,fastmath=False)
def numba_dot2(a, b, c, d):
  ret = 0. + 0j
  for i in range(len(a)):
    ret += a[i] * b[i] + c[i] * d[i]
  return ret

# this is marginally faster than just doing the dots the honest way
@numba.jit(nopython=True,fastmath=False)
def numba_dot3(an,bn,pin,tin):
  ret1 = 0. + 0j
  ret2 = 0. + 0j
  for i in range(len(an)):
    an0 = an[i]
    bn0 = bn[i]
    pin0 = pin[i]
    tin0 = tin[i]
    ret1 += an0 * pin0 + bn0 * tin0
    ret2 += an0 * tin0 + bn0 * pin0
  return (ret1, ret2)

def mie_S12(coeffs,u):
  #return mie_S12old(coeffs,u) # use this to fall back to non-numba version
  return mie_S12_backend(coeffs.nmax, coeffs.an, coeffs.bn, u)

  # use these for self-tests
  #s1, s2 = mie_S12_backend(coeffs.nmax, coeffs.an, coeffs.bn, u)
  #return (np.complex128(s1), np.complex128(s2)) #required to comply with the self-test 

def mie_S12_pt(coeffs,pin, tin):
  return mie_S12_backend_pt(coeffs.nmax, coeffs.an, coeffs.bn, pin, tin)

@numba.jit(nopython=True)
def mie_S12_backend(nmax,an,bn,u):
    """
    Do note that pin and tin do not depend on the refractive index, and thus
    can be shared between runs of different mr, mi
    """
    p = mie_p(u, nmax)
    t = mie_t(u, nmax, p)
    n2 = [float(2 * i + 1) / (i * (i + 1)) for i in range(1,nmax+1)]
    for i in range(nmax):
      p[i] = p[i] * n2[i]
      t[i] = t[i] * n2[i]
    return mie_S12_backend_pt(nmax,an,bn,p,t)

@numba.jit(nopython=True)
def mie_S12_backend_pt(nmax,an,bn,pin, tin):
    """The amplitude scattering matrix.
    """

    # for very large particles this function consumes a lot of time. Possibly the summation is not very efficient?
    """
    this is called extremely many times in most large iterations (once for each angle) and thus
    there is a lot of rationale in optimizing it as much as possible
    """

    # dot2 and dot3 are slightly faster than regular dot, about 5%
    s1 = numba_dot(an,pin)+numba_dot(bn,tin)
    s2 = numba_dot(an,tin)+numba_dot(bn,pin)
    #s1 = numba_dot2(an,pin,bn,tin)
    #s2 = numba_dot2(an,tin,bn,tin)
    return (s1, s2)
    #return numba_dot3(an,bn,pin,tin)

def mie_S12old(coeffs,u):
    """The amplitude scattering matrix.
    """
    (pin,tin) = mie_ptold(u,coeffs.nmax)
    n = arange(1, coeffs.nmax+1, dtype=float)
    n2 = (2*n+1)/(n*(n+1))
    pin *= n2
    tin *= n2

    s1 = dot(coeffs.an,pin)+dot(coeffs.bn,tin)
    s2 = dot(coeffs.an,tin)+dot(coeffs.bn,pin)
    return (s1, s2)

# issues with getting an explicit signature to work, so use the automatic
@numba.jit(nopython=True)
def mie_p(u, nmax):
    p = [0. for i in range(nmax)]
    p[0] = 1
    p[1] = 3*u
    nn = [float(i) for i in range(2,nmax)]

    for n in nn:
        n_i = int(n)
        p[n_i] = (2*n+1)/n*p[n_i-1]*u - (n+1)/n*p[n_i-2]

    return p


"""
Numba version is faster than numpy version
"""
@numba.jit(nopython=True,fastmath=False)
def mie_t(u, nmax, p):
    nn = [float(i) for i in range(2,nmax)]
    t = [0. for i in range(nmax)]
    t[0] = u
    t[1] = 6*u**2 - 3
    for n in nn:
        n_i = int(n)
        t[n_i] = (n+1) * u * p[n_i] - (n+2) * p[n_i-1]
    return t


def mie_pt(u, nmax):
  return mie_ptnumba(u, nmax)

def mie_ptold(u,nmax):
    u = float(u)
    p = zeros(nmax, dtype=float)
    p[0] = 1
    p[1] = 3*u
    t = zeros(nmax, dtype=float)
    t[0] = u
    t[1] = 6*u**2 - 3

    nn = arange(2,nmax,dtype=float)

    for n in nn:
        n_i = int(n)
        p[n_i] = (2*n+1)/n*p[n_i-1]*u - (n+1)/n*p[n_i-2]
    
    t[2:] = (nn+1)*u*p[2:] - (nn+2)*p[1:-1]

    return (p,t)

@numba.jit(nopython=True)
def mie_ptnumba(u,nmax):
    u = float(u)

    p = mie_p(u, nmax)
    t = mie_t(u, nmax, p)

    n2 = [float(2 * i + 1) / (i * (i + 1)) for i in range(1,nmax+1)]
    # faster to just store these n2-multiplied terms, since they only depend on nmax?
    # it is massively faster for size 0.01 ... 1000 range (this function is about 10x faster, total about 33%)
    for i in range(nmax):
      p[i] = p[i] * n2[i]
      t[i] = t[i] * n2[i]

    return (array(p),array(t))
