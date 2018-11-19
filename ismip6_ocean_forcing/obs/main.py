import os
import xarray
import numpy

from ismip6_ocean_forcing.obs import extrap
from ismip6_ocean_forcing.woa.main import process_woa
from ismip6_ocean_forcing.meop.main import process_meop
from ismip6_ocean_forcing.en4.main import process_en4
from ismip6_ocean_forcing.thermal_forcing.main import compute_thermal_forcing
from ismip6_ocean_forcing.remap.res import get_res


def process_obs(config):
    '''
    process and combine WOA, MEOP, EN4 obs and extrapolate into regions with
    missing data
    '''

    if not config.getboolean('observations', 'compute'):
        return

    try:
        os.makedirs('obs')
    except OSError:
        pass

    resFinal = get_res(config, extrap=False)

    startYear = 1995
    endYear = 2017
    decades = '{:04d}-{:04d}'.format(startYear, endYear)

    process_woa(config, decades)
    process_meop(config)
    process_en4(config, startYear, endYear)

    print('Combining and extrapolating observations...')

    _combine_obs(config, decades)

    extrap.extrap_obs(config, decades)

    tempFileName = \
        'obs/obs_temperature_{}_{}.nc'.format(decades, resFinal)
    salinFileName = \
        'obs/obs_salinity_{}_{}.nc'.format(decades, resFinal)
    outFileName = \
        'obs/obs_thermal_forcing_{}_{}.nc'.format(decades, resFinal)
    compute_thermal_forcing(tempFileName, salinFileName, outFileName)

    print('  Done.')


def _combine_obs(config, decades):
    print('  Combining observations...')
    res = get_res(config, extrap=True)
    for fieldName in ['temperature', 'salinity']:
        outFileName = 'obs/obs_{}_{}_{}.nc'.format(fieldName, decades, res)
        if os.path.exists(outFileName):
            continue
        datasets = {}
        datasets['woa'] = xarray.open_dataset('woa/woa_{}_{}_{}.nc'.format(
                    fieldName, decades, res))
        datasets['meop'] = xarray.open_dataset('meop/meop_{}_{}.nc'.format(
                fieldName, res))
        datasets['en4'] = xarray.open_dataset('en4/en4_{}_{}_{}.nc'.format(
                    fieldName, decades, res))
        ds = datasets['en4'].drop(fieldName)
        nx = ds.sizes['x']
        ny = ds.sizes['y']
        nz = ds.sizes['z']
        field = numpy.zeros((nz, ny, nx))
        weights = numpy.zeros(field.shape)
        for name in datasets:
            localField = datasets[name][fieldName].values
            mask = numpy.isfinite(localField)
            field[mask] += localField[mask]
            weights[mask] += 1.

        mask = weights > 0
        field[mask] /= weights[mask]
        field[numpy.logical_not(mask)] = numpy.nan
        ds[fieldName] = (('z', 'y', 'x'), field)
        ds[fieldName].attrs = datasets['en4'][fieldName].attrs

        ds.to_netcdf(outFileName)

    print('    Done.')

