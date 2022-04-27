#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup

long_description = """A Python code for computing the scattering properties
of single- and dual-layered spheres with an easy-to-use object oriented
interface.

Based on code by C. Mätzler; ported and published with permission.

Requires NumPy and SciPy.
"""

setup(name='pymiecoated',
      version='0.2.0',
      download_url=\
          'https://github.com/jleinonen/pymiecoated/archive/v0.2.0.tar.gz',
      description='Single- and dual-layered Mie scattering computations',
      author='Jussi Leinonen',
      author_email='jsleinonen@gmail.com',
      url='http://code.google.com/p/pymiecoated/',
      packages=['pymiecoated','pymiecoated.demos','pymiecoated.test'],
      license='MIT',
      long_description = long_description,
     )
