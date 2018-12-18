#!/usr/bin/env python
import os
import xarray
import numpy

climFirstYear = 1995
climLastYear = 2014
climYears = '{:04d}-{:04d}'.format(climFirstYear, climLastYear)

units = {'temperature': 'Degrees C',
         'salinity': 'PSU',
         'thermal_forcing': 'Degrees C'}

longNames = {'temperature': 'In Situ Seawater Temperature',
             'salinity': 'Seawater Salinity',
             'thermal_forcing': 'Thermal Forcing (in situ Temperature minus '
                                'freezing temperature)'}

for firstYear, lastYear in [(1850, 1994), (1995, 2100)]:
    outFolder = '{:04d}-{:04d}'.format(firstYear, lastYear)
    try:
        os.makedirs(outFolder)
    except OSError:
        pass
    for field in ['temperature', 'salinity', 'thermal_forcing']:
        outFileName = '{}/CCSM4_{}_8km_x_60m.nc'.format(outFolder, field)
        print(outFileName)
        inFileNames = []
        for first in range(firstYear, lastYear, 20):
            last = min(first + 19, lastYear)
            inFileNames.append('{:04d}-{:04d}/anomaly_{:04d}-{:04d}_plus_obs/'
                               'CCSM4_{}_8km_x_60m.nc'.format(
                                   first, last, climFirstYear, climLastYear,
                                   field))
        ds = xarray.open_mfdataset(inFileNames)
        ds[field] = ds[field].astype(numpy.float32)
        ds[field].attrs['units'] = units[field]
        ds[field].attrs['long_name'] = longNames[field]
        ds.to_netcdf(outFileName)
