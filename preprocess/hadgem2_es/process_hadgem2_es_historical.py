#!/usr/bin/env python

import argparse
import xarray
import numpy
import os
import warnings

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', dest='out_dir', metavar='DIR', default='.',
                    type=str, help='output directory')
args = parser.parse_args()


def compute_yearly_mean(inFileNames):
    # crop to below 48 S and take annual mean from monthly data
    print('{}...{} to {}'.format(inFileNames[0], inFileNames[-1], outFileName))

    ds = xarray.open_mfdataset(inFileNames, combine='nested', concat_dim='time')
    dsFirst = xarray.open_dataset(inFileNames[0])

    # crop to Southern Ocean
    ds = ds.isel(lat=slice(0, 43))

    for coord in ['lev_bnds', 'lon_bnds', 'lat_bnds']:
        ds.coords[coord] = dsFirst[coord]

    # annual mean
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ds = ds.groupby('time.year').mean('time', keep_attrs=True)

    # convert back to CF-compliant time
    ds = ds.rename({'year': 'time'})
    ds['time'] = 365.0*ds.time
    ds.time.attrs['bounds'] = "time_bnds"
    ds.time.attrs['units'] = "days since 0000-01-01 00:00:00"
    ds.time.attrs['calendar'] = "noleap"
    ds.time.attrs['axis'] = "T"
    ds.time.attrs['long_name'] = "time"
    ds.time.attrs['standard_name'] = "time"

    timeBounds = numpy.zeros((ds.sizes['time'], 2))
    timeBounds[:, 0] = ds.time.values
    timeBounds[:, 1] = ds.time.values + 365.0
    ds['time_bnds'] = (('time', 'bnds'), timeBounds)

    return ds


dates = ['185912-186911',
         '186912-187911',
         '187912-188911',
         '188912-189911',
         '189912-190911',
         '190912-191911',
         '191912-192911',
         '192912-193911',
         '193912-194911',
         '194912-195911',
         '195912-196911',
         '196912-197911',
         '197912-198911',
         '198912-199911',
         '199912-200512']

fileNames = {}
for field in ['so', 'thetao']:
    fileNames[field] = []
    for date in dates:
        inFileName = '{}/{}_Omon_HadGEM2-ES_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)
        fileNames[field].append(inFileName)

dates = ['200512-201511']

for field in ['so', 'thetao']:
    for date in dates:
        inFileName = '{}/{}_Omon_HadGEM2-ES_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)
        fileNames[field].append(inFileName)


for field in ['so', 'thetao']:
    outFileName = \
        '{}/{}_annual_HadGEM2-ES_rcp85_r1i1p1_186001-201412.nc'.format(
            args.out_dir, field)
    if os.path.exists(outFileName):
        continue

    ds = compute_yearly_mean(fileNames[field])

    # we don't want 1859 or 2015
    ds = ds.isel(time=slice(1, ds.sizes['time']-1))
    ds.to_netcdf(outFileName)
