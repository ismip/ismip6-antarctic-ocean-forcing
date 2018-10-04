import os

from ismip6_ocean_forcing.model import remap
from ismip6_ocean_forcing.extrap.horiz import make_3D_bed_mask, extrap_horiz, \
    extrap_grounded_above_sea_level
from ismip6_ocean_forcing.extrap.vert import extrap_vert
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res
from ismip6_ocean_forcing.remap.interp1d import remap_vertical


def extrapolate_model(config):
    '''
    Extrapolate CMIP model temperature and salinity into ice-shelf cavities and
    IMBIE basins. Ten, compute thermal forcing.
    '''

    if not config.getboolean('model', 'compute'):
        return

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

    resExtrap = get_res(config, extrap=True)
    resFinal = get_res(config, extrap=False)
    hres = get_horiz_res(config)
    modelName = config.get('model', 'name')

    inFileName = '{}/{}_temperature_{}.nc'.format(
        modelFolder, modelName, resExtrap)
    bedMaskFileName = '{}/bed_mask_{}.nc'.format(
        modelFolder, resExtrap)
    bedFileName = 'bedmap2/bedmap2_{}.nc'.format(hres)
    basinNumberFileName = 'imbie/basinNumbers_{}.nc'.format(hres)

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    for fieldName in ['temperature', 'salinity']:
        inFileName = '{}/{}_{}_{}.nc'.format(
            modelFolder, modelName, fieldName, resExtrap)
        outFileName = '{}/{}_{}_{}_extrap_horiz.nc'.format(
            modelFolder, modelName, fieldName, resExtrap)

        progressDir = '{}/progress_{}'.format(modelFolder, fieldName)
        matrixDir = '{}/matrices'.format(modelName.lower())
        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir)

    for fieldName in ['temperature', 'salinity']:
        inFileName = '{}/{}_{}_{}_extrap_horiz.nc'.format(
            modelFolder, modelName, fieldName, resExtrap)

        outFileName = '{}/{}_{}_{}_extrap_vert.nc'.format(
            modelFolder, modelName, fieldName, resExtrap)

        extrap_vert(config, inFileName, outFileName, fieldName)

    inFileNames = {}
    outFileNames = {}
    for fieldName in ['temperature', 'salinity']:

        inFileNames[fieldName] = '{}/{}_{}_{}_extrap_vert.nc'.format(
            modelFolder, modelName, fieldName, resExtrap)
        outFileNames[fieldName] = '{}/{}_{}_{}_extrap_vert.nc'.format(
            modelFolder, modelName, fieldName, resFinal)

    remap_vertical(config, inFileNames, outFileNames, extrap=False)

    for fieldName in ['temperature', 'salinity']:

        inFileName = '{}/{}_{}_{}_extrap_vert.nc'.format(
            modelFolder, modelName, fieldName, resFinal)

        outFileName = '{}/{}_{}_{}.nc'.format(
            modelFolder, modelName, fieldName, resFinal)

        progressDir = '{}/progress_{}'.format(modelFolder, fieldName)
        matrixDir = '{}/matrices'.format(modelName.lower())
        extrap_grounded_above_sea_level(config, inFileName, outFileName,
                                        fieldName, progressDir, matrixDir)
