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

    # make all non-time variables coordinates before computing mean
    for var in ds.data_vars:
        if 'time' not in ds[var].dims:
            ds.coords[var] = ds[var]

    # crop to Southern Ocean
    ds = ds.isel(lat=slice(0, 36))

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

    # add the depth coordinate
    # for k <= nsigma: z = eta + sigma*(min(depth_c,depth)+eta)
    # for k > nsigma: z = zlev
    zIndex = xarray.DataArray.from_dict({'dims': ('lev',),
                                         'data': numpy.arange(ds.sizes['lev'])})
    zIndex.coords['lev'] = ds.lev
    z = ds.eta + ds.sigma*(numpy.minimum(ds.depth_c, ds.depth) + ds.eta)
    z = z.transpose('time', 'lev', 'lat', 'lon')
    z = z.where(zIndex < ds.nsigma, -ds.zlev)

    z_bnds = ds.eta + ds.sigma_bnds*(numpy.minimum(ds.depth_c, ds.depth) \
             + ds.eta)
    z_bnds = z_bnds.transpose('time', 'lev', 'lat', 'lon', 'bnds')
    z_bnds = z_bnds.where(zIndex < ds.nsigma, -ds.zlev_bnds)

    ds['z_full'] = z
    ds.z_full.attrs['units'] = 'm'
    ds.z_full.attrs['description'] = \
        'vertical coordinate at cell centers and mid level'
    ds.z_full.attrs['bounds'] = 'z_bnds'
    ds.z_full.attrs['axis'] = 'Z'

    ds['z_full_bnds'] = z_bnds
    ds.z_full_bnds.attrs['units'] = 'm'
    ds.z_full_bnds.attrs['description'] = \
        'vertical coordinate at cell centers and level interfaces'

    ds = ds.rename({'lev': 'z'})

    ds['z'] = -ds.zlev
    ds.z.attrs['units'] = 'm'
    ds.z.attrs['description'] = \
        'approximate z-level vertical coordinate'
    ds.z.attrs['bounds'] = 'zlev_bnds'
    ds.z.attrs['axis'] = 'Z'

    ds.coords['z_bnds'] = -ds.zlev_bnds

    ds = ds.drop(['eta', 'depth_c', 'depth', 'sigma', 'nsigma', 
                  'sigma_bnds', 'zlev', 'zlev_bnds', 'lev_bnds'])

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
        inFileName = '{}/{}_Omon_MIROC-ESM-CHEM_historical_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_MIROC-ESM-CHEM_historical_r1i1p1_{}.nc'.format(
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
        inFileName = '{}/{}_Omon_MIROC-ESM-CHEM_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        outFileName = '{}/{}_annual_MIROC-ESM-CHEM_rcp85_r1i1p1_{}.nc'.format(
            args.out_dir, field, date)

        compute_yearly_mean(inFileName, outFileName)

for field in ['so', 'thetao']:
    outFileName = '{}/{}_annual_MIROC-ESM-CHEM_rcp85_r1i1p1_185001-210012.nc'.format(
        args.out_dir, field)
    if not os.path.exists(outFileName):
        print(outFileName)

        # combine it all into a single data set
        ds = xarray.open_mfdataset('{}/{}_annual_MIROC-ESM-CHEM_*_r1i1p1_*.nc'.format(
            args.out_dir, field), concat_dim='time')
        ds.to_netcdf(outFileName)
