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
    dsIn = dsIn.rename({'z_t': 'lev', 'd2': 'bnds', 'TLONG': 'lon',
                        'TLAT': 'lat'})
    ds = xarray.Dataset()
    if 'SALT' in dsIn:
        ds['so'] = dsIn['SALT']
    if 'TEMP' in dsIn:
        ds['thetao'] = dsIn['TEMP']
    ds.coords['lev'] *= 0.01
    ds.lev.attrs['units'] = 'm'
    ds.lev.attrs['bounds'] = 'lev_bnds'
    lev_bnds = numpy.zeros((dsIn.sizes['lev'], dsIn.sizes['bnds']))
    lev_bnds[:, 0] = 0.01*dsIn.z_w_top
    lev_bnds[:, 1] = 0.01*dsIn.z_w_bot
    ds.coords['lev_bnds'] = (('lev', 'bnds'), lev_bnds)
    ds = ds.drop(['ULONG', 'ULAT'])

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
run = 'e21_f09_g17'

dates = {'thetao': ['195001-199912',
                    '200001-201412'],
         'so': ['195001-199912',
                '200001-201412']}

inField = {'thetao': 'TEMP', 'so': 'SALT'}

histFiles = {}
for field in dates:
    histFiles[field] = []
    for date in dates[field]:
        inFileName = '{}/b.e21.BHIST.f09_g17.CMIP6-historical.011.pop.h.' \
                     '{}.{}.nc'.format(args.out_dir, inField[field], date)

        tempFileName = '{}/temp.nc'.format(args.out_dir)
        outFileName = '{}/{}_annual_{}_historical_{}_{}.nc'.format(
            args.out_dir, field, model, run, date)

        compute_yearly_mean(inFileName, outFileName, tempFileName)
        histFiles[field].append(outFileName)

dates = {'thetao': ['201501-206412',
                    '206501-210012'],
         'so': ['201501-206412',
                '206501-210012']}
for scenario in ['ssp585']:
    scenarioFiles = {}
    for field in dates:
        scenarioFiles[field] = []
        for date in dates[field]:
            inFileName = '{}/b.e21.BSSP585cmip6.f09_g17.CMIP6-SSP5-8.5.002.' \
                         'pop.h.{}.{}.nc'.format(args.out_dir, inField[field],
                                                 date)

            tempFileName = '{}/temp.nc'.format(args.out_dir)

            outFileName = '{}/{}_annual_{}_{}_{}_{}.nc'.format(
                args.out_dir, field, model, scenario, run, date)

            compute_yearly_mean(inFileName, outFileName, tempFileName)
            scenarioFiles[field].append(outFileName)

    for field in ['so', 'thetao']:
        outFileName = \
            '{}/{}_annual_{}_{}_{}_199501-210012.nc'.format(
                args.out_dir, field, model, scenario, run)
        if not os.path.exists(outFileName):
            print(outFileName)

            # combine it all into a single data set
            ds = xarray.open_mfdataset(histFiles[field] + scenarioFiles[field],
                                       combine='nested', concat_dim='time',
                                       use_cftime=True)

            mask = ds['time.year'] >= 1995
            tIndices = numpy.nonzero(mask.values)[0]
            ds = ds.isel(time=tIndices)
            encoding = {'time': {'units': 'days since 0000-01-01'}}
            ds.to_netcdf(outFileName, encoding=encoding)
