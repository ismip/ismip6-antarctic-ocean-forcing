#!/usr/bin/env python

import argparse
import xarray
import numpy
import os
import warnings

parser = argparse.ArgumentParser(
    description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-o', dest='out_dir', metavar='DIR', required=True,
                    type=str, help='output directory')
args = parser.parse_args()


def compute_yearly_mean(inFileName, outFileName, tempFileName):
    # crop to below 48 S and take annual mean over the data
    if os.path.exists(outFileName):
        return
    print('{} to {}'.format(inFileName, outFileName))

    dsIn = xarray.open_dataset(inFileName, chunks={'time': 24})
    dsIn = dsIn.rename({'d2': 'bnds'})
    ds = xarray.Dataset()
    if 'so' in dsIn:
        ds['so'] = dsIn['so']
    if 'thetao' in dsIn:
        ds['thetao'] = dsIn['thetao']
    ds.coords['lev'] *= 0.01
    ds.lev.attrs['units'] = 'm'
    ds.lev.attrs['bounds'] = 'lev_bnds'
    ds.coords['lev_bnds'] = dsIn.lev_bnds

    # crop to Southern Ocean
    minLat = ds.lat.min(dim='nlon')
    mask = minLat <= -48.
    yIndices = numpy.nonzero(mask.values)[0]
    ds = ds.isel(nlat=yIndices)

    # write out and read back to get a clean start
    ds.to_netcdf(tempFileName)
    ds.close()
    ds = xarray.open_dataset(tempFileName)

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
    timeBounds[:, 0] = ds.time.values
    timeBounds[:, 1] = ds.time.values + 365
    ds['time_bnds'] = (('time', 'bnds'), timeBounds)

    ds = xarray.decode_cf(ds, use_cftime=True)

    encoding = {'time': {'units': 'days since 0000-01-01'}}
    ds.to_netcdf(outFileName, encoding=encoding)


model = 'CESM2'
run = 'r11i1p1f1'

dates = {'thetao': ['185001-189912',
                    '190001-194912',
                    '195001-199912',
                    '200001-201412'],
         'so': ['185001-189912',
                '190001-194912',
                '195001-199912',
                '200001-201412']}

histFiles = {}
for field in dates:
    histFiles[field] = []
    for date in dates[field]:
        inFileName = '{}/{}_Omon_{}_historical_{}_gn_{}.nc'.format(
            args.out_dir, field, model, run, date)

        tempFileName = '{}/temp.nc'.format(args.out_dir)
        outFileName = '{}/{}_annual_{}_historical_{}_{}.nc'.format(
            args.out_dir, field, model, run, date)

        compute_yearly_mean(inFileName, outFileName, tempFileName)
        histFiles[field].append(outFileName)

run = 'r2i1p1f1'
dates = {'thetao': [],
         'so': []}
for scenario in ['ssp585']:
    scenarioFiles = {}
    for field in dates:
        scenarioFiles[field] = []
        for date in dates[field]:
            inFileName = '{}/{}_Omon_{}_{}_{}_gn_{}.nc'.format(
                args.out_dir, field, model, scenario, run, date)

            tempFileName = '{}/temp.nc'.format(args.out_dir)

            outFileName = '{}/{}_annual_{}_{}_{}_{}.nc'.format(
                args.out_dir, field, model, scenario, run, date)

            compute_yearly_mean(inFileName, outFileName, tempFileName)
            scenarioFiles[field].append(outFileName)

run = 'r11i1p1f1'
for scenario in ['historical']:
    for field in ['so', 'thetao']:
        outFileName = \
            '{}/{}_annual_{}_{}_{}_185001-201412.nc'.format(
                args.out_dir, field, model, scenario, run)
        if not os.path.exists(outFileName):
            print(outFileName)

            # combine it all into a single data set
            ds = xarray.open_mfdataset(histFiles[field] + scenarioFiles[field],
                                       combine='nested', concat_dim='time',
                                       use_cftime=True)

            mask = ds['time.year'] <= 2014
            tIndices = numpy.nonzero(mask.values)[0]
            ds = ds.isel(time=tIndices)
            encoding = {'time': {'units': 'days since 0000-01-01'}}
            ds.to_netcdf(outFileName, encoding=encoding)
