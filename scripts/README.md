Table of Contents
=================

* [GEOSmie](#geosmie)
   * [v2.0.0](#v2.0.0)

<!-- Created by https://github.com/ekalinin/github-markdown-toc -->


# GEOSmie

See [Kemppinen et al. 2022](https://gmao.gsfc.nasa.gov/pubs/docs/Kemppinen1447.pdf) for 
documentation of the approach. In the following is some brief document of the particlular
characteristics of specific optical tables generated, sorted by a carefully thought
out versioning scheme!

## V2.0.0

This version is intended to be as close as possible to the legacy IDL-based tables, which 
in turn are mainly based on OPAC-like assumptions of particle properties (growth factor
and complex refractive index). Some details are made explicit in Chin et al. 2002 and 
Colarco et al. 2010. Specific details are provided in the JSON files in the geosparticles
sub-directory and additional information in Kemppinen et al. 2022.

### Carbonaceous Aerosols
Following Chin et al. 2002 the particle size distributions are lognormal modes specified
in terms of number mode median radius (in microns) with a given width and truncated at
0.3 micron radius (dry particle). Growth factors specify the ratio of wet-to-dry particle
radius at the given RH. Refractive indices are given spectrally for the dry particles
and are volume weighted with water refractive index to give an effective refractive index
for the Mie calculations. For each of the carbonaceous species there are two size modes:
the first mode is hydrophobic (undergoes no growth with RH) and the second mode is
hydrophilic (grows with RH with specified growth factor).

| Species | rnum [um] | sigma | Density [kg m-3] | Notes    |
| ---     | ---       | ---   | ---              | ---      |
| BC      | 0.0118    | 2.0   | 1000             | bc.json  |
| OC      | 0.0212    | 2.20  | 1800             | oc.json  |
| BR      | 0.0212    | 2.20  | 1800             | brc.json |

### Sulfate
Following Chin et al. 2002 the particle size distributions are lognormal modes specified
in terms of number mode median radius (in microns) with a given width and truncated at
0.3 micron radius (dry particle). Growth factors specify the ratio of wet-to-dry particle
radius at the given RH. Refractive indices are given spectrally for the dry particles
and are volume weighted with water refractive index to give an effective refractive index
for the Mie calculations. A single mode of sulfate is provided here.

| Species | rnum [um] | sigma | Density [kg m-3] | Notes    |
| ---     | ---       | ---   | ---              | ---      |
| SU      | 0.0695    | 2.03  | 1000             | su.json  |


