#!/usr/bin/env python

from setuptools import setup, find_packages

install_requires = ['dask',
                    'descartes>=1.1.0',
                    'gsw',
                    'imageio',
                    'matplotlib',
                    'netcdf4',
                    'numpy',
                    'requests',
                    'progressbar2',
                    'pyshp',
                    'scikit-fmm',
                    'scipy',
                    'shapely>=2.0',
                    'xarray']

version = '1.0'

setup(name='ismip6_ocean_forcing',
      version=version,
      description='A python package for creating ISMIP6 ocean forcing data '
                  'sets.',
      url='https://github.com/xylar/ismip6-ocean-forcing',
      author='Xylar Asay-Davis',
      author_email='xylarstorm@gmail.com',
      license='MIT',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'License :: OSI Approved :: MIT License',
          'Operating System :: OS Independent',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering',
      ],
      packages=find_packages(),
      package_data={'ismip6_ocean_forcing': ['default.cfg']},
      install_requires=install_requires,
      entry_points={'console_scripts':
                    ['ismip6_ocean_forcing = ismip6_ocean_forcing.__main__:main']})
