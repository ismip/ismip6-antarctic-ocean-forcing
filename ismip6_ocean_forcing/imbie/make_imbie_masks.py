import os
import zipfile

from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.imbie import images, extend
from ismip6_ocean_forcing.remap.res import get_horiz_res


def make_imbie_masks(config):
    '''
    Download geojson files defining IMBIE basins and create masks that extend
    these basins into the open ocean
    '''
    res = get_horiz_res(config)

    try:
        os.makedirs('imbie')
    except OSError:
        pass

    if os.path.exists('imbie/basinNumbers_{}.nc'.format(res)):
        # we're already done
        return

    print('Building IMBIE basin masks...')

    basins = {'A-Ap': ['A-Ap'],
              'Ap-B': ['Ap-B'],
              'B-C': ['B-C'],
              'C-Cp': ['C-Cp'],
              'Cp-D': ['Cp-D'],
              'D-Dp': ['D-Dp'],
              'Dp-E': ['Dp-E'],
              'E-F': ['E-Ep', 'Ep-F'],
              'F-G': ['F-G'],
              'G-H': ['G-H'],
              'H-Hp': ['H-Hp'],
              'Hp-I': ['Hp-I'],
              'I-Ipp': ['I-Ipp'],
              'Ipp-J': ['Ipp-J'],
              'J-K': ['J-Jpp', 'Jpp-K'],
              'K-A': ['K-A']}

    bedFileName = 'bedmap2/bedmap2_{}.nc'.format(res)

    _download_imbie()

    images.write_basin_images(res, bedFileName, basins)

    extend.extend_imbie_masks(res, basins, bedFileName)
    print('  Done.')


def _download_imbie():
    '''
    Download the geojson files that define the IMBIE basins
    '''

    if not os.path.exists(
            'imbie/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'):
        # download
        geojsonURL = 'http://imbie.org/wp-content/uploads/2016/09/'

        fileNames = ['ANT_Basins_IMBIE2_v1.6.zip']

        download_files(fileNames, geojsonURL, 'imbie')

        print('Decompressing IMBIE2 data...')
        # unzip
        with zipfile.ZipFile('imbie/ANT_Basins_IMBIE2_v1.6.zip', 'r') as f:
            f.extractall('imbie/ANT_Basins_IMBIE2_v1.6/')
        print('  Done.')
