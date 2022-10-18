# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed

- Fixed a bug where one-line kernel shape distribution files were read incorrectly

### Added

- Add code and runscript (convertkernels.py and runkernelconversion.py) for converting GRASP-like kernel format for a NetCDF file that GEOSmie can use as kernels. Also included are two preset (in data/kernelconversion/) parameter files for reading two kernel types used so far.
- Read kernel shape distribution from a JSON-supplied file
- Add JSON files for several experimental dust types, included in geosparticles/experimental/ directory
- Add band wavelength upper/lower edge information to band-averaged files

### Changed

- Kernel file path moved to a kerneal\_params dict, breaking backward compatibility with old non-Mie JSON files
- Changed the default behavior of generate\_geos\_optics.sh to not run dust, given dust needs extra steps to run

### Removed

### Deprecated

## [1.0.0] - 2022-09-28

### Added

- Initialized repository with usual markdown files (README, CHANGELOG, etc.)
- Add changelog enforcer
- Add required script and data files for the package
- Add running instructions in README
