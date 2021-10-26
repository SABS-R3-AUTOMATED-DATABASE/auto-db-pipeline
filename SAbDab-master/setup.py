#!/usr/bin/env python
"""
setup.py
Setup script for ABDB.
"""
from distutils.core import setup, Extension

# Specify numpy include dir for compiling KDTree
#import numpy as np
#numpy_include_dir = np.get_include()

setup(name     ="ABDB",
        version='1.1',
        description='',
        author='James Dunbar, Konrad Krawczyk, Jinwoo Leem',
        author_email='leem@stats.ox.ac.uk',
        packages=[
          'ABDB',
          'ABDB.AB_Utils',
          'ABDB.Annotate',
          'ABDB.AbPDB',
          'ABDB.ABDB_updater',
          'ABDB.ABangle',

      ],
        package_dir={
          'ABDB':'lib/python/ABDB',
          'ABDB.AB_Utils': 'lib/python/ABDB/AB_Utils',
          'ABDB.Annotate': 'lib/python/ABDB/Annotate',
          'ABDB.AbPDB': 'lib/python/ABDB/AbPDB',
          'ABDB.ABDB_updater': 'lib/python/ABDB/ABDB_updater',
          'ABDB.ABangle': 'lib/python/ABDB/ABangle',

      },
        package_data = 
            {
                'ABDB.ABangle': ['dat/*']
            },
        scripts = ['bin/SAbDab']
     )

