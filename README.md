Table of Contents
=================

* [GEOSmie](#geosmie)
   * [How to build GEOSmie](#how-to-build-geosmie)
   * [Main code](#main-code)
      * [Calculations at individual wavelengths](#calculations-at-individual-wavelengths)
      * [Calculations over wavelength bands](#calculations-over-wavelength-bands)
   * [GSF expansion code](#gsf-expansion-code)
   * [Kernel conversion code](#kernel-conversion-code)

<!-- Created by https://github.com/ekalinin/github-markdown-toc -->


# GEOSmie

The GEOSmie package is used to generate aerosol optical property lookup tables (LUTs)
for use with the GOCART2G package in the GEOS Earth system model and related tools. It
supports table generation either from direct Mie calculations or from integration over
externally provided kernels in a specified format.

In general, two sorts of LUTs are produced for each of the aerosol species simulated:
- The first is the monochromatic table (i.e., optical properties at specified wavelengths)
which has a form like `optics_XX.vY.Y.Y.nc4` where `XX` is a label indicating the 
aerosol species and `vY.Y.Y` is some sort of versioning string (e.g., `optics_DU.v2.0.0.nc4`).
This file is produced with a call to `runoptics.py` and a subsequent call to `rungsf.py`
to add the expansion moments of the phase function elements. The monochromatic tables are
used for higher precision calculation of aerosol optical properties at specified
wavelengths, for example the aerosol optical depth at 550 nm.
- The second sort of LUT is generated from the monochromatic table and is averaged over
spectral bands dictated by the GEOS-internal radiative transfer code used for radiative
flux and forcing calculation (i.e., RRTMG). The tables have a form like `optics_XX.vY.Y.Y.RRTMG.nc4`
(e.g., `optics_DU.v2.0.0.RRTMG.nc4`). The wavelength index in that table is related
to the radiative transfer code band indices in that case. These files are produced by a
call to `runbands.py`. See below for usage examples.

See [Kemppinen et al. 2022](https://gmao.gsfc.nasa.gov/pubs/docs/Kemppinen1447.pdf) for 
documentation of the approach. Below is a brief summary of the structure and usage of the
package. More specific information about the optical property assumptions is provided in
the `README.md` in the `scripts` subdirectory, which also contains version-specific 
table generation scripts in an attempt to make this as turnkey as possible.

The package consists of the following parts:

- Mie code (pymiecoated/)
- Main code (root directory), which can calculate single-scattering and bulk optical properties for both individual wavelengths and wavelength bands
- Generalized spherical function expansion code (gsf/)
- Kernel generation code (root directory)

## How to Build GEOSmie
### Preliminary Steps

#### Load Build Modules

In your `.bashrc` or `.tcshrc` or other rc file add a line:

##### NCCS

```
module use -a /discover/swdev/gmao_SIteam/modulefiles-SLES15
```

##### NAS
```
module use -a /nobackup/gmao_SIteam/modulefiles
```

##### GMAO Desktops

On the GMAO desktops, the SI Team modulefiles should automatically be
part of running `module avail` but if not, they are in:

```
module use -a /ford1/share/gmao_SIteam/modulefiles
```

Also do this in any interactive window you have. This allows you to get module files needed to correctly checkout and build the model.

Now load the `GEOSenv` module:
```
module load GEOSenv
```
which obtains the latest `git`, `CMake`, etc. modules needed to build.

##### Non-GMAO Desktops

Target is for the 614 lab cluster, assumed to have `git`, `cmake`, and `mepo` installed.
It is not necessary (or useful, or--really--possible) in this case to fuss about with 
loading the `GEOSenv` module. The only notable modification needed then is to check
`cmake_it` code and change the fortran reference as necessary; e.g., change `ifort` to
`gfortran` on bender.

### Use mepo to clone the repository

[Mepo](https://github.com/GEOS-ESM/mepo) is a multiple repository tool available on github.

```
mepo clone git@github.com:GEOS-ESM/GEOSmie.git
```

##### Slow clones

If you notice your clone is taking a while, we recommend running:

```
mepo config set clone.partial blobless
```

This is a one-time command that tells mepo to use blobless clones for all future clones. 
Blobless clones are much faster than the default clone method, 
especially for repositories with a large history like MAPL.

#### Build the Model

##### Step 1: Load Compiler, MPI Stack, and Baselibs

On tcsh:
```
source env@/g5_modules
```
or on bash:
```
source env@/g5_modules.sh
```

##### Create Build Directory

##### Step 2: Run CMake

CMake generates the Makefiles needed to build the model.

```
./cmake_it
```

This will build the code in `build/` and will install in `install`.

NOTE: You can choose any directory for your build and install in `cmake_it`. That
said, what you pass to `--install-prefix=` *MUST* be a full path.

###### Building with Debugging Flags

To build with debugging flags add:

```
-DCMAKE_BUILD_TYPE=Debug
```
to the cmake line.

##### Step 3: Build and Install with CMake

```
cmake --build build --target install -j6
```

If you put your build in a directory other than `build` use that in the
above command.

## Setup an experiment

Navigate to `install/bin`. Run the `geosmie_setup` script and follow the prompts:

```
./geosmie_setup.py
```

This will create an independent experiment directory for running GEOSmie.
You can choose between starting scripts found in `src/scripts` as a template.
These scripts implement the following ways to run GEOSmie.

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

To run the conversion code:

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
./runkernelconversion.py --filename data/kernelconversion/filename.json --dest /path/to/output/directory
```

where filename.json is a parameter file that specifies various aspects on how the kernels should be read, and has to be updated if the incoming kernel format changes in virtually any way. Sample JSON files are provided for GRASP kernels (grasp.json) and GRASP-like Saito et al. kernels (saito.json).

## Contributing

Please check out our [contributing guidelines](CONTRIBUTING.md).

## License

All files are currently licensed under the Apache-2.0 license, see [`LICENSE`](LICENSE).

