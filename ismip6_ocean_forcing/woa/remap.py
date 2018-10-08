import xarray
import os
import numpy

from ismip6_ocean_forcing.remap.interp1d import weights_and_indices, \
    interp_depth

from ismip6_ocean_forcing.remap.descriptor import get_antarctic_descriptor, \
    get_lat_lon_descriptor

from ismip6_ocean_forcing.remap.remapper import Remapper
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res


def remap_woa(config):
    _combine_climatologies()
    _interp_z(config)
    _remap(config)


def _combine_climatologies():

    bothExist = True
    for field in ['t', 's']:
        outFileName = 'woa/woa13_95B2_{}00_04v2_no_time.nc'.format(field)
        if not os.path.exists(outFileName):
            bothExist = False
            break
    if bothExist:
        return

    print('  Combining WOA climatologies...')
    for field in ['t', 's']:
        outFileName = 'woa/woa13_95B2_{}00_04v2_no_time.nc'.format(field)
        if os.path.exists(outFileName):
            continue
        ds = xarray.Dataset()
        for decade, weight in [('95A4', 10./18.), ('A5B2', 8./18.)]:
            fileName = 'woa/woa13_{}_{}00_04v2.nc'.format(decade, field)
            print('    {}'.format(fileName))
            dsIn = xarray.open_dataset(fileName, decode_times=False)
            depth = dsIn['depth']
            dsIn['depth'] = -depth
            dsIn['depth'].attrs = depth.attrs
            dsIn['depth'].attrs['positive'] = 'up'
            depth = dsIn['depth_bnds']
            dsIn['depth_bnds'] = -depth
            dsIn['depth_bnds'].attrs = depth.attrs
            for coord in ['lon', 'lat', 'depth']:
                if coord not in ds.coords:
                    ds[coord] = dsIn[coord]
                bounds = '{}_bnds'.format(coord)
                if bounds not in ds.data_vars:
                    ds[bounds] = dsIn[bounds]
            varName = '{}_mn'.format(field)
            contribution = weight*dsIn[varName].isel(time=0)
            if varName in ds.data_vars:
                ds[varName] += contribution
            else:
                ds[varName] = contribution
                ds[varName].attrs = dsIn[varName].attrs

        ds = ds.drop('time')
        ds.to_netcdf(outFileName)


def _interp_z(config):

    bothExist = True
    for field in ['t', 's']:
        outFileName = 'woa/woa13_95B2_{}00_04v2_interp_z.nc'.format(field)
        if not os.path.exists(outFileName):
            bothExist = False
            break
    if bothExist:
        return

    print('  Interpolate in depth to common grid...')
    dz = config.getfloat('grid', 'dzExtrap')
    nz = config.getint('grid', 'nzExtrap')
    zOut = dz*numpy.arange(nz+1)

    for field in ['t', 's']:
        inFileName = 'woa/woa13_95B2_{}00_04v2_no_time.nc'.format(field)
        outFileName = 'woa/woa13_95B2_{}00_04v2_interp_z.nc'.format(field)
        print('    {}'.format(outFileName))
        dsIn = xarray.open_dataset(inFileName)
        zIn = numpy.zeros(dsIn.sizes['depth']+1)
        zIn[0:-1] = dsIn.depth_bnds[:, 0]
        zIn[-1] = dsIn.depth_bnds[-1, 1]
        weights, inIndices = weights_and_indices(xInBounds=zIn,
                                                 xOutBounds=zOut,
                                                 xDim='z')

        varName = '{}_mn'.format(field)
        dsOut = xarray.Dataset()
        result = interp_depth(dsIn[varName], weights, inIndices,
                              normalizationThreshold=0.1)
        dsOut[varName] = result
        for attrName in ['units', 'standard_name', 'long_name']:
            dsOut[varName].attrs[attrName] = dsIn[varName].attrs[attrName]
        for coord in ['lon', 'lat']:
            dsOut[coord] = dsIn[coord]
            bounds = '{}_bnds'.format(coord)
            dsOut[bounds] = dsIn[bounds]
        z = 0.5*(zOut[0:-1] + zOut[1:])
        z_bnds = numpy.zeros((len(z), 2))
        z_bnds[:, 0] = zOut[0:-1]
        z_bnds[:, 1] = zOut[1:]
        dsOut['z'] = (('z',), z)
        dsOut.z.attrs = dsIn.depth.attrs
        dsOut.z.attrs['bounds'] = 'z_bnds'
        dsOut['z_bnds'] = (('z', 'nbounds'), z_bnds)
        dsOut.z_bnds.attrs = dsIn.depth_bnds.attrs
        dsOut[varName].coords['z'] = dsOut.z

        dsOut.to_netcdf(outFileName)


def _remap(config):

    res = get_res(config)
    hres = get_horiz_res(config)
    bothExist = True
    for fieldName in ['temperature', 'salinity']:
        outFileName = 'woa/woa_{}_1995-2012_{}.nc'.format(fieldName, res)
        if not os.path.exists(outFileName):
            bothExist = False
            break
    if bothExist:
        return

    print('  Remapping to {} grid...'.format(res))
    for field, fieldName in [['t', 'temperature'], ['s', 'salinity']]:
        inFileName = 'woa/woa13_95B2_{}00_04v2_interp_z.nc'.format(field)
        outGridFileName = 'ismip6/{}_grid.nc'.format(hres)
        outFileName = 'woa/woa_{}_1995-2012_{}.nc'.format(fieldName, res)
        print('    {}'.format(outFileName))

        varName = '{}_mn'.format(field)

        inDescriptor = get_lat_lon_descriptor(inFileName)
        outDescriptor = get_antarctic_descriptor(outGridFileName)

        mappingFileName = 'woa/map_{}_to_{}.nc'.format(inDescriptor.meshName,
                                                       outDescriptor.meshName)

        remapper = Remapper(inDescriptor, outDescriptor, mappingFileName)

        remapper.build_mapping_file(method='bilinear')

        ds = xarray.open_dataset(inFileName)
        ds = ds.rename({varName: fieldName})

        dsOut = remapper.remap(ds, renormalizationThreshold=0.1)

        for attrName in ['units', 'standard_name', 'long_name']:
            dsOut[fieldName].attrs[attrName] = ds[fieldName].attrs[attrName]
        dsOut.z.attrs = ds.z.attrs

        dsOut.to_netcdf(outFileName)
