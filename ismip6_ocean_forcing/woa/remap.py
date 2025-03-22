import xarray
import os
import numpy
import gsw

from pyremap import get_polar_descriptor_from_file, LatLonGridDescriptor, \
    Remapper

from ismip6_ocean_forcing.remap.interp1d import weights_and_indices, \
    interp_depth
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res


def remap_woa(config, woaDecades, woaWeights, decades):
    _combine_climatologies(woaDecades, woaWeights, decades)
    _interp_z(config, decades)
    _remap(config, decades)


def _combine_climatologies(woaDecades, woaWeights, decades):

    bothExist = True
    for fieldName in ['temperature', 'salinity']:
        outFileName = 'woa/woa18_{}_{}_no_time.nc'.format(decades, fieldName)
        if not os.path.exists(outFileName):
            bothExist = False
            break
    if bothExist:
        return

    print('  Combining WOA climatologies...')
    for fieldName, shortName, inVarName in [('temperature', 't', 't_mn'),
                                            ('salinity', 's', 's_mn')]:
        outFileName = 'woa/woa18_{}_{}_no_time.nc'.format(decades, fieldName)
        if os.path.exists(outFileName):
            continue
        ds = xarray.Dataset()
        weightSum = None
        outField = None
        for decade, weight in zip(woaDecades, woaWeights):
            fileName = 'woa/woa18_{}_{}00_04.nc'.format(decade, shortName)
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
            contribution = weight*dsIn[inVarName].isel(time=0)
            mask = contribution.notnull().values
            if weightSum is None:
                weightSum = weight*mask
                outField = numpy.zeros(mask.shape)
                outField[mask] = contribution.values[mask]
            else:
                weightSum += weight*mask
                outField[mask] += contribution.values[mask]

        mask = numpy.logical_and(
            numpy.isfinite(weightSum),
            weightSum > 0
        )
        outField[mask] /= weightSum[mask]
        outField[numpy.logical_not(mask)] = numpy.nan
        ds[fieldName] = (('depth', 'lat', 'lon'), outField)
        ds[fieldName].attrs = dsIn[inVarName].attrs

        ds.to_netcdf(outFileName)


def _interp_z(config, decades):

    bothExist = True
    for fieldName in ['temperature', 'salinity']:
        outFileName = 'woa/woa18_{}_{}_interp_z.nc'.format(decades, fieldName)
        if not os.path.exists(outFileName):
            bothExist = False
            break
    if bothExist:
        return

    print('  Interpolate in depth to common grid...')
    dz = config.getfloat('grid', 'dzExtrap')
    nz = config.getint('grid', 'nzExtrap')
    zOut = dz*numpy.arange(nz+1)

    for fieldName in ['temperature', 'salinity']:
        inFileName = 'woa/woa18_{}_{}_no_time.nc'.format(decades, fieldName)
        outFileName = 'woa/woa18_{}_{}_interp_z.nc'.format(decades, fieldName)
        print('    {}'.format(outFileName))
        dsIn = xarray.open_dataset(inFileName)
        zIn = numpy.zeros(dsIn.sizes['depth']+1)
        zIn[0:-1] = dsIn.depth_bnds[:, 0]
        zIn[-1] = dsIn.depth_bnds[-1, 1]
        weights, inIndices = weights_and_indices(xInBounds=zIn,
                                                 xOutBounds=zOut,
                                                 xDim='z')

        varName = fieldName
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


def _remap(config, decades):

    res = get_res(config)
    hres = get_horiz_res(config)
    bothExist = True
    for fieldName in ['temperature', 'salinity']:
        outFileName = 'woa/woa_{}_{}_{}.nc'.format(fieldName, decades, res)
        if not os.path.exists(outFileName):
            bothExist = False
            break
    if bothExist:
        return

    print('  Remapping to {} grid...'.format(res))
    for fieldName in ['temperature', 'salinity']:
        inFileName = 'woa/woa18_{}_{}_interp_z.nc'.format(decades, fieldName)
        outGridFileName = 'ismip6/{}_grid.nc'.format(hres)
        outFileName = 'woa/woa_{}_{}_{}.nc'.format(fieldName, decades, res)
        if os.path.exists(outFileName):
            continue
        print('    {}'.format(outFileName))

        varName = fieldName

        inDescriptor = LatLonGridDescriptor.read(fileName=inFileName)
        outDescriptor = get_polar_descriptor_from_file(outGridFileName,
                                                       projection='antarctic')

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


def _in_situ_to_potential_temperature(dsTemp, dsSalin):
    z = dsTemp.z.values
    lat = dsTemp.lat.values
    lon = dsTemp.lon.values

    nz = len(z)
    ny, nx = lat.shape

    dsPotTemp = dsTemp.drop('in_situ_temperature')
    pt = numpy.nan*numpy.ones((nz, ny, nx))
    for zIndex in range(nz):
        pressure = gsw.p_from_z(z[zIndex], lat)
        in_situ_temp = dsTemp.in_situ_temperature[zIndex, :, :].values
        salin = dsSalin.salinity[zIndex, :, :].values
        mask = numpy.isfinite(in_situ_temp)
        SA = gsw.SA_from_SP(salin[mask], pressure[mask], lon[mask],
                            lat[mask])
        ptSlice = pt[zIndex, :, :]
        ptSlice[mask] = gsw.pt_from_t(SA, in_situ_temp[mask], pressure[mask],
                                      p_ref=0.)
        pt[zIndex, :, :] = ptSlice

    dsPotTemp['temperature'] = (('z', 'y', 'x'), pt)
    dsPotTemp['temperature'].attrs = dsTemp.in_situ_temperature.attrs

    return dsPotTemp

