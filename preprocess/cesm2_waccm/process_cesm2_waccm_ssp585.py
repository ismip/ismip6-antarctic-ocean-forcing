#!/usr/bin/env python

import argparse
import xarray
import numpy
import os
import warnings

from dask.diagnostics import ProgressBar


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

    with xarray.open_dataset(inFileName, chunks={'time': 24}) as dsIn:
        lev_bnds = dsIn.lev_bnds.values
        minLat = dsIn.lat.min(dim='nlon')
        mask = minLat <= -48.
        dsIn = dsIn.where(mask, drop=True)
        # write out and read back to get a clean start
        delayed_obj = dsIn.to_netcdf(tempFileName, compute=False)
        with ProgressBar():
            results = delayed_obj.compute()

    dsIn = xarray.open_dataset(tempFileName)

    dsIn = dsIn.rename({'d2': 'bnds'})
    ds = xarray.Dataset()
    if 'so' in dsIn:
        ds['so'] = dsIn['so']
    if 'thetao' in dsIn:
        ds['thetao'] = dsIn['thetao']
    ds.coords['lev'] = 0.01*dsIn.coords['lev']
    ds.lev.attrs['units'] = 'm'
    ds.lev.attrs['bounds'] = 'lev_bnds'
    ds.coords['lev_bnds'] = (('lev', 'bnds'), lev_bnds)
    ds.lev_bnds.attrs['units'] = 'm'

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
    delayed_obj = ds.to_netcdf(outFileName, encoding=encoding, compute=False)
    with ProgressBar():
        results = delayed_obj.compute()


model = 'CESM2-WACCM'
run = 'r1i1p1f1'

dates = {'thetao': ['185001-201412'],
         'so': ['185001-201412']}

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

run = 'r1i1p1f1'
dates = {'thetao': ['201501-210012',
                    '210101-215012',
                    '215101-220012',
                    '220101-225012',
                    '225101-229912'],
         'so': ['201501-210012',
                '210101-215012',
                '215101-220012',
                '220101-225012',
                '225101-229912']}

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

    for field in ['so', 'thetao']:
        outFileName = \
            '{}/{}_annual_{}_{}_{}_185001-229912.nc'.format(
                args.out_dir, field, model, scenario, run)
        if not os.path.exists(outFileName):
            print(outFileName)

            # combine it all into a single data set
            ds = xarray.open_mfdataset(histFiles[field] + scenarioFiles[field],
                                       combine='nested', concat_dim='time',
                                       use_cftime=True)

            encoding = {'time': {'units': 'days since 0000-01-01'}}
            ds.to_netcdf(outFileName, encoding=encoding)
