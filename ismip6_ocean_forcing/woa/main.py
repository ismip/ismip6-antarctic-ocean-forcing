import os
from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.woa import remap, extrap


def extrapolate_woa(config):
    '''
    Download WOA 13 v2 temperature and salinity, and extrapolate them into
    ice-shelf cavities and IMBIE basins
    '''
    try:
        os.makedirs('woa')
    except OSError:
        pass

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

    print('  Done.')
