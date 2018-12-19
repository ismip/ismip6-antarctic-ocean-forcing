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


dates = ['185001-185312',
         '185401-185712',
         '185801-186112',
         '186201-186512',
         '186601-186912',
         '187001-187312',
         '187401-187712',
         '187801-188112',
         '188201-188512',
         '188601-188912',
         '189001-189312',
         '189401-189712',
         '189801-190112',
         '190201-190512',
         '190601-190912',
         '191001-191312',
         '191401-191712',
         '191801-192112',
         '192201-192512',
         '192601-192912',
         '193001-193312',
         '193401-193712',
         '193801-194112',
         '194201-194512',
         '194601-194912',
         '195001-195312',
         '195401-195712',
         '195801-196112',
         '196201-196512',
         '196601-196912',
         '197001-197312',
         '197401-197712',
         '197801-198112',
         '198201-198512',
         '198601-198912',
         '199001-199312',
         '199401-199712',
         '199801-200112',
         '200201-200512']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_NorESM1-M_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_NorESM1-M_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName)

dates = ['200601-200912',
         '201001-201312',  
         '201401-201712',  
         '201801-202112',  
         '202201-202512',  
         '202601-202912',  
         '203001-203312',  
         '203401-203712',  
         '203801-204112',  
         '204201-204512',  
         '204601-204912',  
         '205001-205312',  
         '205401-205712',  
         '205801-206112',  
         '206201-206512',  
         '206601-206912',  
         '207001-207312',  
         '207401-207712',  
         '207801-208112',  
         '208201-208512',  
         '208601-208912',  
         '209001-209312',  
         '209401-209712',  
         '209801-210012']

for date in dates:
    for field in ['so', 'thetao']:
        inFileName = '{}/{}_Omon_NorESM1-M_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_NorESM1-M_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName)

for field in ['so', 'thetao']:
    outFileName = '{}/{}_annual_NorESM1-M_rcp85_r1i1p1_185001-210012.nc'.format(
        args.out_dir, field)
    if not os.path.exists(outFileName):
        print(outFileName)

        # combine it all into a single data set
        ds = xarray.open_mfdataset('{}/{}_annual_NorESM1-M_*_r1i1p1_*.nc'.format(
            args.out_dir, field), concat_dim='time')
        ds.to_netcdf(outFileName)
