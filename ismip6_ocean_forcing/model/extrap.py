import os

from ismip6_ocean_forcing.model import remap
from ismip6_ocean_forcing.extrap.horiz import make_3D_bed_mask, extrap_horiz
from ismip6_ocean_forcing.extrap.vert import extrap_vert
from ismip6_ocean_forcing.remap.res import get_res


def extrapolate_model(config):
    '''
    Extrapolate CMIP model temperature and salinity into ice-shelf cavities and
    IMBIE basins. Ten, compute thermal forcing.
    '''

    modelName = config.get('model', 'name')
    subfolder = config.get('model', 'folder')
    modelFolder = '{}/{}'.format(modelName.lower(), subfolder)

    try:
        os.makedirs(modelFolder)
    except OSError:
        pass

    print('Extrapolate {} model results...'.format(modelName))

    remap.remap_model(config, modelFolder)

    _extrap_model(config, modelFolder)

    print('  Done.')


def _extrap_model(config, modelFolder):

    res = get_res(config)
    modelName = config.get('model', 'name')

    inFileName = '{}/{}_temperature_{}.nc'.format(
        modelFolder, modelName, res)
    bedMaskFileName = '{}/bed_mask_{}.nc'.format(
        modelFolder, res)
    bedFileName = 'bedmap2/bedmap2_{}.nc'.format(res)
    basinNumberFileName = 'imbie/basinNumbers.nc'

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    for fieldName in ['temperature', 'salinity']:
        inFileName = '{}/{}_{}_{}.nc'.format(
            modelFolder, modelName, fieldName, res)
        outFileName = '{}/{}_{}_{}_extrap_horiz.nc'.format(
            modelFolder, modelName, fieldName, res)

        progressDir = '{}/progress_{}'.format(modelFolder, fieldName)
        matrixDir = '{}/matrices'.format(modelName.lower())
        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir)

    for fieldName in ['temperature', 'salinity']:
        inFileName = '{}/{}_{}_{}_extrap_horiz.nc'.format(
            modelFolder, modelName, fieldName, res)

        outFileName = '{}/{}_{}_{}_extrap_vert.nc'.format(
            modelFolder, modelName, fieldName, res)

        extrap_vert(config, inFileName, outFileName, fieldName)
