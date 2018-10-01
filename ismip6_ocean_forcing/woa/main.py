import os
from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.woa import remap, extrap
from ismip6_ocean_forcing.thermal_forcing.main import compute_thermal_forcing
from ismip6_ocean_forcing.remap.res import get_res


def extrapolate_woa(config):
    '''
    Download WOA 13 v2 temperature and salinity, and extrapolate them into
    ice-shelf cavities and IMBIE basins
    '''
    try:
        os.makedirs('woa')
    except OSError:
        pass

    resFinal = get_res(config, extrap=False)

    print('Extrapolate World Ocean Atlas...')

    baseURL = 'https://data.nodc.noaa.gov/thredds/fileServer/woa/WOA13/' \
              'DATAv2/temperature/netcdf/95A4/0.25/'
    fileNames = ['woa13_95A4_t00_04v2.nc']
    download_files(fileNames, baseURL, 'woa')
    baseURL = 'https://data.nodc.noaa.gov/thredds/fileServer/woa/WOA13/' \
              'DATAv2/salinity/netcdf/95A4/0.25/'
    fileNames = ['woa13_95A4_s00_04v2.nc']
    download_files(fileNames, baseURL, 'woa')
    baseURL = 'https://data.nodc.noaa.gov/thredds/fileServer/woa/WOA13/' \
              'DATAv2/temperature/netcdf/A5B2/0.25/'
    fileNames = ['woa13_A5B2_t00_04v2.nc']
    download_files(fileNames, baseURL, 'woa')
    baseURL = 'https://data.nodc.noaa.gov/thredds/fileServer/woa/WOA13/' \
              'DATAv2/salinity/netcdf/A5B2/0.25/'
    fileNames = ['woa13_A5B2_s00_04v2.nc']
    download_files(fileNames, baseURL, 'woa')

    remap.remap_woa(config)

    extrap.extrap_woa(config)

    tempFileName = \
        'woa/woa_temperature_1995-2012_{}.nc'.format(resFinal)
    salinFileName = \
        'woa/woa_salinity_1995-2012_{}.nc'.format(resFinal)
    outFileName = \
        'woa/woa_thermal_forcing_1995-2012_{}.nc'.format(resFinal)
    compute_thermal_forcing(tempFileName, salinFileName, outFileName)

    print('  Done.')
