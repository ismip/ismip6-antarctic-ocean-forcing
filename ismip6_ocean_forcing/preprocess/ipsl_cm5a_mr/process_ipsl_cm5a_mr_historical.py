#!/usr/bin/env python

import argparse
import xarray
import numpy
import os
import warnings
import glob


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

    for coord in ['lev_bnds', 'lon_vertices', 'lat_vertices']:
        ds.coords[coord] = ds[coord]

    # crop to Southern Ocean
    ds = ds.where(ds.lat <= -48., drop=True)

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

    timeBounds = numpy.zeros((ds.sizes['time'], ds.sizes['bnds']))
    timeBounds[:, 0] = 365.0*ds.time.values
    timeBounds[:, 1] = 365.0*(ds.time.values+1)
    ds['time_bnds'] = (('time', 'bnds'), timeBounds)

    ds.to_netcdf(outFileName)


dates = ['185001-189912',
         '190001-194912',
         '195001-199912',
         '200001-200512']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_IPSL-CM5A-MR_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = \
            '{}/{}_annual_IPSL-CM5A-MR_historical_r1i1p1_{}.nc'.format(
                args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName)

dates = ['200601-205512']

for scenario in ['rcp26', 'rcp85']:
    for date in dates:
        for field in ['so', 'thetao']:
            inFileName = '{}/{}_Omon_IPSL-CM5A-MR_{}_r1i1p1_{}.nc'.format(
                args.out_dir, field, scenario, date)

            outFileName = '{}/{}_annual_IPSL-CM5A-MR_{}_r1i1p1_{}.nc'.format(
                args.out_dir, field, scenario, date)

            compute_yearly_mean(inFileName, outFileName)

for scenario in ['rcp26', 'rcp85']:
    for field in ['so', 'thetao']:
        outFileName = \
            '{}/{}_annual_IPSL-CM5A-MR_{}_r1i1p1_185001-201412.nc'.format(
                args.out_dir, field, scenario)
        if not os.path.exists(outFileName):
            print(outFileName)

            # combine it all into a single data set
            files = list()
            for expt in ['historical', scenario]:
                files = files + glob.glob(
                    '{}/{}_annual_IPSL-CM5A-MR_{}_r1i1p1_*.nc'.format(
                         args.out_dir, field, expt))

            ds = xarray.open_mfdataset(files, combine='nested',
                                       concat_dim='time', use_cftime=True)
            mask = ds['time.year'] <= 2014
            tIndices = numpy.nonzero(mask.values)[0]
            ds = ds.isel(time=tIndices)
            encoding = {'time': {'units': 'days since 0000-01-01'}}
            ds.to_netcdf(outFileName)
