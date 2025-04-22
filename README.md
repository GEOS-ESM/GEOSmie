Table of Contents
=================

* [GEOSmie](#geosmie)
   * [Mie code](#mie-code)
   * [Main code](#main-code)
      * [Calculations at individual wavelengths](#calculations-at-individual-wavelengths)
      * [Calculations over wavelength bands](#calculations-over-wavelength-bands)
   * [GSF expansion code](#gsf-expansion-code)
   * [Kernel conversion code](#kernel-conversion-code)

<!-- Created by https://github.com/ekalinin/github-markdown-toc -->


# GEOSmie

The GEOSmie package is used to generate aerosol optical property lookup tables (LUTs)
for use with the GOCART2G package in the GEOS Earth system model and related tools.

In general, two sorts of LUTs are produced for each of the aerosol species simulated.
The first is the monochromatic table (i.e., optical properties at specified wavelengths)
which has a form like "optics_XX.vY.Y.Y.nc4" where "XX" is a label indicating the 
aerosol species and "vY.Y.Y" is some sort of versioning string (e.g., optics_DU.v2.2.2.nc4).
This file is produced with a call to "runoptics.py" and a subsequent call to "rungsf.py"
to add the expansion moments of the phase function elements. The monochromatic tables are
used for higher precision calculation of aerosol optical properties at specified
wavelengths, for example the aerosol optical depth at 550 nm.
The second sort of LUT is generated from the monochromatic table and is averaged over
spectral bands dictated by the GEOS-internal radiative transfer code used for radiative
flux and forcing calculation (i.e., RRTMG). The wavelength index in that table is related
to the radiative transfer code band indices in that case. These files are produced by a
call to "runbands.py". See below for usage examples.

See [Kemppinen et al. 2022](https://gmao.gsfc.nasa.gov/pubs/docs/Kemppinen1447.pdf) for 
documentation of the approach. Below is a brief summary of the structure and usage of the
package. More specific information about the optical property assumptions is provided in
the README.md in the `scripts` subdirectory, which also contains version-specific 
table generation scripts in an attempt to make this as turnkey as possible.

The package consists of the following parts:

- Mie code (pymiecoated/)
- Main code (root directory), which can calculate single-scattering and bulk optical properties for both individual wavelengths and wavelength bands
- Generalized spherical function expansion code (gsf/)

## Mie code

Before any spherical aerosol calculations can be done the Mie code needs to be installed.

Starting in the root directory of this repository:

```bash
cd pymiecoated
python setup.py install
```

In some systems the following may be needed instead:

```bash
python setup.py install --user
```

The Mie code does not need to be used directly; runoptics.py calls it as needed. The functions within the module can be used directly to run custom Mie simulations but that is beyond the scope of this README.


## Main code

The main single-scattering property calculator consists of two parts: calculation at individual wavelengths, and calculation over wavelength bands. The first part is always necessary, the second needs to be run only if wavelength band integrated files are necessary.

### Calculations at individual wavelengths

To calculate single-scattering properties at individual wavelengths:

```bash
./runoptics.py --name path/to/file.json
```

Where file.json defines all of the parameters of the calculation. ".json" can be omitted from the runoptics.py command as a convenience. See JSON examples under geosparticles/ for current GEOS aerosol particles. For users who only want single-scattering properties at individual wavelengths this is all that is needed, and further processing is not needed.

runoptics.py has a few additional options, use

```bash
./runoptics.py --help
```

to see them

### Calculations over wavelength bands

To calculate single-scattering properties at wavelength bands:

```bash
./runbands.py --filename filename.nc
```

where filename.nc should be an output file from runoptics.py

For additional options, see

```bash
./runbands.py --help
```

## GSF expansion code

The generalized spherical function expansion code is needed to convert the phase functions in the optical tables to be compatible with the VLIDORT radiative transfer software. It is not needed unless you are running VLIDORT directly.

The code uses spher_expan.f code by Michael Mishchenko. Before conversions can be done the Fortran code needs to be compiled. 

Starting in the root directory of this repository:

```bash
cd gsf
gfortran spher_expan.f
```

Other compilers beyond gfortran may work. Testing has been performed with gfortran 11.2.0 on Gentoo.

To run the conversion code, assuming the code has been compiled and you are in the gsf/ subdirectory:

```bash
./rungsf.py --filename ../filename.nc
```

where ../filename.nc should be an output file from runoptics.py

For additional options, see

```bash
./rungsf.py --help
```

## Kernel conversion code

External kernels (e.g. GRASP spheroids) have to be converted into a specific NetCDF format before GEOSmie can use them. For details on the format, see the full documentation.

A tool is provided for converting kernels from GRASP-like format to the GEOSmie-compatible NetCDF file.

```bash
python runkernelconversion.py --filename data/kernelconversion/filename.json --dest /path/to/output/directory
```

where filename.json is a parameter file that specifies various aspects on how the kernels should be read, and has to be updated if the incoming kernel format changes in virtually any way. Sample JSON files are provided for GRASP kernels (grasp.json) and GRASP-like Saito et al. kernels (saito.json).
