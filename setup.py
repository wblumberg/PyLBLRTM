#!/usr/bin/env python

from distutils.core import setup

setup(name='PyLBLRTM',
      version='1.0',
      description='Python scripts to run and analyze the LBLRTM output',
      author='Greg Blumberg',
      author_email='wblumberg@ou.edu',
      packages=['pylblrtm'],
      package_data = {"": ["*.npy", "*.md", "*.txt"]}

)


