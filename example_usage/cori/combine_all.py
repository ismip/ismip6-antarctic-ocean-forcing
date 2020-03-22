#!/usr/bin/env python
import os
import xarray
import numpy

def combine(start, end, model, climFirstYear=1995, climLastYear=2014):

    climYears = '{:04d}-{:04d}'.format(climFirstYear, climLastYear)
    
    units = {'temperature': 'Degrees C',
             'salinity': 'PSU',
             'thermal_forcing': 'Degrees C'}
    
    longNames = {'temperature': 'In Situ Seawater Temperature',
                 'salinity': 'Seawater Salinity',
                 'thermal_forcing': 'Thermal Forcing (in situ Temperature minus'
                                    ' freezing temperature)'}
    
    for firstYear, lastYear in [(start, climFirstYear-1), (climFirstYear, end)]:
        outFolder = '{:04d}-{:04d}'.format(firstYear, lastYear)
        try:
            os.makedirs(outFolder)
        except OSError:
            pass
        for field in ['temperature', 'salinity', 'thermal_forcing']:
            outFileName = '{}/{}_{}_8km_x_60m.nc'.format(outFolder, model,
                                                         field)
            print(outFileName)
            inFileNames = []
            for first in range(firstYear, lastYear, 20):
                last = min(first + 19, lastYear)
                inFileNames.append(
                   '{:04d}-{:04d}/anomaly_{:04d}-{:04d}_plus_obs/'
                   '{}_{}_8km_x_60m.nc'.format(first, last, climFirstYear,
                                               climLastYear, model, field))
            ds = xarray.open_mfdataset(inFileNames, combine='nested',
                                       concat_dim='time')
            ds[field] = ds[field].astype(numpy.float32)
            ds[field].attrs['units'] = units[field]
            ds[field].attrs['long_name'] = longNames[field]
            ds.to_netcdf(outFileName)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", dest="start", type=int,
                    help="start year of the time series")
    parser.add_argument("-e", dest="end", type=int,
                    help="end year of the time series")
    parser.add_argument("-m", dest="model", type=str,
                    help="name of the model")
    args = parser.parse_args() 

    setup_years(args.start, args.end, args.model, args.scenario, args.ensemble)
    setup_decades(args.start, args.end, args.model)

