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


def compute_yearly_mean(inFileName, outFileName):
    # crop to below 48 S and take annual mean from monthly data
    if os.path.exists(outFileName):
        return
    print('{} to {}'.format(inFileName, outFileName))

    ds = xarray.open_dataset(inFileName)

    # crop to Southern Ocean
    ds = ds.isel(lat=slice(0, 44))

    for coord in ['lev_bnds', 'lon_bnds', 'lat_bnds']:
        ds.coords[coord] = ds[coord]

    # annual mean
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        ds = ds.groupby('time.year').mean('time', keep_attrs=True)

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
    timeBounds[:, 0] = 365.0*ds.time.values
    timeBounds[:, 1] = 365.0*(ds.time.values+1)
    ds['time_bnds'] = (('time', 'bnds'), timeBounds)

    ds.to_netcdf(outFileName)


dates = ['199001-199912',
         '200001-200512']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_CSIRO-Mk3-6-0_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_CSIRO-Mk3-6-0_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName)

dates = ['200601-201512',
         '201601-202512',
         '202601-203512',
         '203601-204512',
         '204601-205512',
         '205601-206512',
         '206601-207512',
         '207601-208512',
         '208601-209512',
         '209601-210012']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_CSIRO-Mk3-6-0_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_CSIRO-Mk3-6-0_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName)

for field in ['so', 'thetao']:
    outFileName = \
        '{}/{}_annual_CSIRO-Mk3-6-0_rcp85_r1i1p1_199001-210012.nc'.format(
            args.out_dir, field)
    if not os.path.exists(outFileName):
        print(outFileName)

        # combine it all into a single data set
        ds = xarray.open_mfdataset(
            '{}/{}_annual_CSIRO-Mk3-6-0_*_r1i1p1_*.nc'.format(
                args.out_dir, field), concat_dim='time')
        ds.to_netcdf(outFileName)
