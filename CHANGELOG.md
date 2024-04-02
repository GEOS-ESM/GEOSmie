# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

### Added

### Changed

### Removed

### Deprecated

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
