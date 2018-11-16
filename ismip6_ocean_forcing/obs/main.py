import os
from ismip6_ocean_forcing.obs import extrap
from ismip6_ocean_forcing.woa.main import process_woa
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

    decades = '1995-2017'

    process_woa(config, decades)

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
    print('    Done.')

