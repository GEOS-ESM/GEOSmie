# GEOSmie
GEOSmie

See [here (TBD)](tbd) for a full documentation. Below is a brief summary of the structure and usage of the package.

The package consists of the following parts:

- Main code (root directory), which can calculate single-scattering and bulk optical properties for both individual wavelengths and wavelength bands
- Mie code (pymiecoated/)
- Generalized spherical function expansion code (gsf/)

## Main code

The main single-scattering property calculator consists of two parts: calculation at individual wavelengths, and calculation over wavelength bands. The first part is always necessary, the second needs to be run only if wavelength band integrated files are necessary.

### Individual wavelengths

To calculate single-scattering properties at individual wavelengths:

python runoptics.py --name path/to/file.json

Where file.json defines all of the parameters of the calculation. ".json" can be omitted from the runoptics.py command as a convenience. See JSON examples under geosparticles/ for current GEOS aerosol particles. runoptics.py has a few additional options, use

python runoptics.py --help

to see them

### Wavelength bands

To calculate single-scattering properties at wavelength bands:

python runbands.py --filename filename.nc

where filename.nc should be an output file from runoptics.py

For additional options, see

python runbands.py --help

## Mie code

Before any spherical aerosol calculations can be done the Mie code needs to be installed.

Starting in the root directory of this repository:

cd pymiecoated
python setup.py install

In some systems the following may be needed instead:

python setup.py install --user

The Mie code does not need to be used directly; runoptics.py calls it as needed. The functions within the module can be used directly but that is beyond the scope of this README.

## GSF expansion code

The generalized spherical function expansion code is needed to convert the phase functions in the optical tables to be compatible with the VLIDORT radiative transfer software. It is not needed unless you are running VLIDORT directly.

The code uses spher_expan.f code by Michael Mishchenko. Before conversions can be done the Fortran code needs to be compiled. 

Starting in the root directory of this repository:

cd gsf
gfortran spher_expan.f

Other compilers beyond gfortran may work. Testing has been performed with gfortran 11.2.0 on Gentoo.

To run the conversion code, assuming the code has been compiled and you are in the gsf/ subdirectory:

python rungsf.py --filename ../filename.nc

where ../filename.nc should be an output file from runoptics.py

For additional options, see

python rungsf.py --help
