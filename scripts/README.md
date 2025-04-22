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

### Carbonaceous Aerosols - BC = Black Carbon, OC = Organic Carbon, BR = Brown Carbon
Following Chin et al. 2002 the particle size distributions are lognormal modes specified
in terms of number mode median radius (in microns) with a given width and truncated at
0.3 micron radius (dry particle). Growth factors specify the ratio of wet-to-dry particle
radius at the given RH. Refractive indices are given spectrally for the dry particles
and are volume weighted with water refractive index to give an effective refractive index
for the Mie calculations. For each of the carbonaceous species there are two size modes:
the first mode is hydrophobic (undergoes no growth with RH) and the second mode is
hydrophilic (grows with RH with specified growth factor). Dry particle refractive indices
are from OPAC (soot00 for BC; waso00 for OC; waso00 is used for BR from wavelengths >= 550 nm
and follow Colarco et al. 2017 for shorter wavelengths).

| Species | rnum [um] | sigma | Density [kg m-3] | Notes    |
| ---     | ---       | ---   | ---              | ---      |
| BC      | 0.0118    | 2.0   | 1000             | bc.json  |
| OC      | 0.0212    | 2.20  | 1800             | oc.json  |
| BR      | 0.0212    | 2.20  | 1800             | brc.json |

### Sulfate - SU
Following Chin et al. 2002 the particle size distributions are lognormal modes specified
in terms of number mode median radius (in microns) with a given width and truncated at
0.3 micron radius (dry particle). Growth factors specify the ratio of wet-to-dry particle
radius at the given RH. Refractive indices are given spectrally for the dry particles using OPAC
(suso00) and are volume weighted with water refractive index to give an effective refractive index
for the Mie calculations. A single mode of sulfate is provided here (su.json).

| Species | rnum [um] | sigma | Density [kg m-3] | Notes    |
| ---     | ---       | ---   | ---              | ---      |
| SU      | 0.0695    | 2.03  | 1000             | su.json  |

### Sea salt - SS
Following Colarco et al. 2010, sea salt aerosol optical properties are computed for each
of five size bins. The v2.0.0 implementation treats the bins as truncated bins with the
optical properties integrated across the bin. Growth factors are based on the Gerber [1986]
parameterization. Properties for the dry particle bins are given below and dry sea salt density 
is assumed to be 2200 kg m-3. A sub-bin distribution is assumed across the bin following Gong [2003].
Refractive indices are given spectrally for the dry particles using OPAC
(sscm00) and are volume weighted with water refractive index to give an effective refractive index
for the Mie calculations.

| Bin | rMin [um] | rMax [um] | Notes   |
| --- | ---       | ---       | ---     |
| 1   | 0.03      | 0.1       | ss.json |
| 2   | 0.1       | 0.5       | ss.json |
| 3   | 0.5       | 1.5       | ss.json |
| 4   | 1.5       | 5.0       | ss.json |
| 5   | 5.0       | 10.0      | ss.json |

### Nitrate - NI
Following Bian et al. 2017, nitrate aerosol optical properties are computed for each of three
size bins (ni.json). The bins are lognormal modes with a specified number mode radius (in microns)
and width (sigma) and are not truncated with sharp edges.
Hygroscopic growth factor is taken to be 1.06 * growth factor of sulfate,
following Fitzgerald 1975. Refractive indices are given spectrally for the dry particles after
Bian et al. 2017 and are volume weighted with water refractive index to give an effective 
refractive index for the Mie calculations.

| Bin | rNum [um] | sigma | density [kg m-3] | Notes   |
| --- | ---       | ---   | ---              | ---     |
| 1   | 0.04455   | 2.03  | 1725             | ni.json |
| 2   | 0.6000    | 2.03  | 2200             | ni.json |
| 3   | 2.3316    | 2.0   | 2650             | ni.json |

### Dust - DU
Dust is novel in that non-spherical optics are considered and Mie calculations are not used.
Optical properties are from the GRASP kernels assuming the spheroidal shape distribution of
Dubovik et al. 2006. Refractive indices are from Colarco et al. 2014, which are a combination
of observation-based shortwave values and OPAC values in the longwave. Five size bins are
considered, but properties are integrated across generally wider distributions than the bin 
edges.

| Bin | rMin [um] | rMax [um] | rNum [uM] | sigma | density [kg m-3] | Notes                            |
| --- | ---       | ---       | ---       | ---   | ---              | ---                              |
| 1   | 0.1       | 1.0       | 0.5173    | 1.327 | 2500             | du-grasp_spheroid-lognormal.json |
| 2   | 1.0       | 1.8       | 1.031     | 1.284 | 2650             | du-grasp_spheroid-lognormal.json |
| 3   | 1.8       | 3.0       | 1.611     | 1.411 | 2650             | du-grasp_spheroid-lognormal.json |
| 4   | 3.0       | 6.0       | 1.198     | 2.028 | 2650             | du-grasp_spheroid-lognormal.json |
| 5   | 6.0       | 10.0      | 2.047     | 2.067 | 2650             | du-grasp_spheroid-lognormal.json |

### Band Average Properties
A simple linear weighting of monochromatic properties is applied across the RRTMG bands (see `bandaverage.py`).

Shortwave bands [cm-1]
| Index | lower | upper |
| ---   | ---   | ---   |
| 1     | 2600  | 3250  |
| 2     | 3250  | 4000  |
| 3     | 4000  | 4650  |
| 4     | 4650  | 5150  |
| 5     | 5150  | 6150  |
| 6     | 6150  | 7700  |
| 7     | 7700  | 8050  |
| 8     | 8050  | 12850 |
| 9     | 12850 | 16000 |
| 10    | 16000 | 22650 |
| 11    | 22650 | 29000 |
| 12    | 29000 | 38000 |
| 13    | 38000 | 50000 |
| 14    | 820   | 2600  |

Longwave bands [cm-1]
| Index | lower | upper |
| ---   | ---   | ---   |
| 15    | 10    | 250   |
| 16    | 250   | 500   |
| 17    | 500   | 630   |
| 18    | 630   | 700   |
| 19    | 700   | 820   |
| 20    | 820   | 980   |
| 21    | 980   | 1080  |
| 22    | 1080  | 1180  |
| 23    | 1180  | 1390  |
| 24    | 1390  | 1480  |
| 25    | 1480  | 1800  |
| 26    | 1800  | 2080  |
| 27    | 2080  | 2250  |
| 28    | 2250  | 2380  |
| 29    | 2380  | 2600  |
| 30    | 2600  | 3250  |

### References
- Bian et al. 2017: https://doi.org/10.5194/acp-17-12911-2017
- Chin et al. 2002: https://journals.ametsoc.org/view/journals/atsc/59/3/1520-0469_2002_059_0461_taotft_2.0.co_2.xml
- Colarco et al. 2010: https://doi.org/10.1029/2009jd012820
- Colarco et al. 2014: https://doi.org/10.1002/2013jd020046
- Colarco et al. 2017: https://doi.org/10.5194/amt-10-4121-2017
- Dubovik et al. 2006: https://doi.org/10.1029/2005jd006619
- Fitzgerald 1975: https://journals.ametsoc.org/view/journals/apme/14/6/1520-0450_1975_014_1044_afftes_2_0_co_2.xml
- Gerber, H. E.: Relative-Humidity Parameterization of the Navy Aerosol Model (NAM), NRL Report 8956, 1 16, 1985.
- OPAC: Hess et al. 1998: https://journals.ametsoc.org/view/journals/bams/79/5/1520-0477_1998_079_0831_opoaac_2_0_co_2.xml
