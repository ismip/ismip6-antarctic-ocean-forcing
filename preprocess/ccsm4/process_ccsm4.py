#!/usr/bin/env python

import argparse
import xarray
import numpy
import os
import warnings

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', dest='out_dir', metavar='DIR',
                    type=str, help='output directory')
args = parser.parse_args()


def compute_yearly_mean(inFileName, outFileName, correctSalinity):
    # crop to below 48 S and take annual mean over the RCP 8.5 data
    if os.path.exists(outFileName):
        return
    print('{} to {}'.format(inFileName, outFileName))

    ds = xarray.open_dataset(inFileName)

    # crop to Southern Ocean
    ds = ds.isel(j=slice(0, 60))

    for coord in ['lev_bnds', 'lon_vertices', 'lat_vertices']:
        ds.coords[coord] = ds[coord]

    # annual mean
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ds = ds.groupby('time.year').mean('time', keep_attrs=True)

    if 'time_bnds' in ds:
        ds = ds.drop('time_bnds')

    # convert back to CF-compliant time
    ds = ds.rename({'year': 'time'})
    ds['time'] = 365.0*ds.time
    ds.time.attrs['bounds'] = "time_bnds"
    ds.time.attrs['units'] = "days since 0000-01-01 00:00:00"
    ds.time.attrs['calendar'] = "noleap"
    ds.time.attrs['axis'] = "T"
    ds.time.attrs['long_name'] = "time"
    ds.time.attrs['standard_name'] = "time"

    timeBounds = numpy.zeros((ds.sizes['time'], ds.sizes['bnds']))
    timeBounds[:, 0] = ds.time.values
    timeBounds[:, 1] = ds.time.values + 365.0
    ds['time_bnds'] = (('time', 'bnds'), timeBounds)

    if field == 'so' and correctSalinity:
        # convert salinity to PSU
        attrs = ds[field].attrs
        ds[field] = 1000.*ds[field]
        ds[field].attrs = attrs
        ds[field].attrs['units'] = 'PSU'

    ds.to_netcdf(outFileName)


dates = ['185001-185912',
         '186001-186912',
         '187001-187912',
         '188001-188912',
         '189001-189912',
         '190001-190912',
         '191001-191912',
         '192001-192912',
         '193001-193912',
         '194001-194912',
         '195001-195912',
         '196001-196912',
         '197001-197912',
         '198001-198912',
         '199001-199912',
         '200001-200512']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_CCSM4_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_CCSM4_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName, correctSalinity=True)

dates = ['200601-200912',
         '201001-201912',
         '202001-202912',
         '203001-203912',
         '204001-204912',
         '205001-205912',
         '206001-206912',
         '207001-207912',
         '208001-208912',
         '209001-210012']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_CCSM4_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_CCSM4_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName, correctSalinity=False)

dates = ['210101-210912',
         '211001-211912',
         '212001-212912',
         '213001-213912',
         '214001-214912',
         '215001-215912',
         '216001-216912',
         '217001-217912',
         '218001-218912',
         '219001-219912',
         '220001-220912',
         '221001-221912',
         '222001-222912',
         '223001-223912',
         '224001-224912',
         '225001-225912',
         '226001-226912',
         '227001-227912',
         '228001-228912',
         '229001-229912',
         '230001-230012']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_CCSM4_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_CCSM4_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName, correctSalinity=True)

for field in ['so', 'thetao']:
    outFileName = '{}/{}_annual_CCSM4_rcp85_r1i1p1_185001-230012.nc'.format(
        args.out_dir, field)
    if not os.path.exists(outFileName):
        print(outFileName)

        # combine it all into a single data set
        ds = xarray.open_mfdataset('{}/{}_annual_CCSM4_*_r1i1p1_*.nc'.format(
            args.out_dir, field), concat_dim='time', decode_times=False)
        ds.to_netcdf(outFileName)
