# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
- fixed eval_gsfun script for legendre overflows and bad indexing of scattering phase functions when plotting
### Added

### Changed

- Updated `.gitignore` to ignore more binary file types as well as the recommended Python patterns

### Removed

### Deprecated

## [1.2.1] - 2025-04-10

### Added

- Added components.yaml and CMakeLists.txt and directory structure to bring
  GMAOpyobs into GEOSmie

## [1.2.0] - 2025-04-10

### Changed

- Changed file output names to end in ".nc4" for consistency with v1.0.0 files
- Changed naming style of band tables for consistency with legacy
- JSON settings files for kernelconversion.py are no longer required to include a "numang"
  key but now do require "renormUpperX", "angFwdPeak", and "kernelScaleFact"
- Normalization of P11 prior to integral calculating g is obtained from GRASP kernel 
  instead of renormalization in kernelconversion.py
- Integral to calculate g uses scipy.integrate to permit variable scattering angle spacing

### Fixed

- Fixed mass extinction efficiency calculation for kernel driven files
- Fixed sign of s12 output
- Scattering matrix elements written by convertkernels.py are now normalized such that P11
  integrates to 4*pi instead of using esoteric GRASP scaling scheme
- Implemented new scheme relying on normalization of P11 in GRASP kernels to prevent
  unphysical decrease in g at large x caused by steep forward scattering peak
- Moved writing of NetCDF inside a with statement so file is not left open on crash
- Ensured nonzero scattering when it's smaller than number of digits in GRASP kernel files

### Added

- Created utils directory with eval_gsfun.py script
- Added scripting directory for version control of table generation
- Long names and descriptions to variables in NetCDF file output by convertkernels.py
- Calculation of graspfactor so kernels with different size parameter spacing can be used 
- Reading of scattering angles from GRASP kernels so 1 deg spacing is not always assumed
- New JSON example files for other GRASP kernel types, including the Saito database 
- Comments at top of convertkernels.py defining keys required in input JSON files
- .gitignore so that .pyc and .nc4 files are not tracked by git


## [1.1.2] = 2024-04-02

### Fixed
- Fixed a bug in the dimension ordering that was only encountered
  for "trivial" RH handling applicable to dust

### Added
- Added a "dust-mie.json" case for spherical dust

## [1.1.1] - 2024-03-21

### Fixed
- Clean up some of the merging of working code with develop

## [1.1.0] - 2024-03-21

### Changed

- All driver files (runoptics, runbands, rungsf) are now "executable"
  from cmd line
- Added some comments to code
- Reordered and renamed dimensions per discussion with modeling group
- Kernel file path moved to a kerneal\_params dict, breaking backward compatibility with old non-Mie JSON files
- Changed the default behavior of generate\_geos\_optics.sh to not run dust, given dust needs extra steps to run

### Fixed

- There were some errors in initial implementation of code to do with
  normalization of hydrated properties that were fixed.
- Fixed a bug where one-line kernel shape distribution files were read incorrectly

### Added

- Add code and runscript (convertkernels.py and runkernelconversion.py) for converting GRASP-like kernel format for a NetCDF file that GEOSmie can use as kernels. Also included are two preset (in data/kernelconversion/) parameter files for reading two kernel types used so far.
- Read kernel shape distribution from a JSON-supplied file
- Add JSON files for several experimental dust types, included in geosparticles/experimental/ directory
- Add band wavelength upper/lower edge information to band-averaged files

## [1.0.1] = 2024-03-21

### Fixed

- Bug fix on integration of backscatter efficiency calculation
  to be consistent with IDL code
- Bug fix on dry mass calculation

## [1.0.0] - 2022-09-28

### Added

- Initialized repository with usual markdown files (README, CHANGELOG, etc.)
- Add changelog enforcer
- Add required script and data files for the package
- Add running instructions in README
