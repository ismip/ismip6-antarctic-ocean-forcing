import os
from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.woa import remap


def process_woa(config, decades):
    '''
    Download WOA 18 temperature and salinity, and remap them to the ISMIP6
    grid
    '''

    try:
        os.makedirs('woa')
    except OSError:
        pass

    print('  Processing World Ocean Atlas...')

    woaDecades = ['95A4', 'A5B7']
    woaWeights = [10., 13.]

    for decade in woaDecades:
        for fieldName, shortName in [('temperature', 't'), ('salinity', 's')]:

            baseURL = 'https://data.nodc.noaa.gov/thredds/fileServer/ncei/' \
                      'woa/{}/{}/0.25/'.format(fieldName, decade)
            fileNames = ['woa18_{}_{}00_04.nc'.format(decade, shortName)]
            download_files(fileNames, baseURL, 'woa')

    remap.remap_woa(config, woaDecades, woaWeights, decades)

    print('    Done.')
