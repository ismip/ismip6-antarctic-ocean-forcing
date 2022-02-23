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
    modelFolder = os.path.join(modelName.lower(), subfolder)

    try:
        os.makedirs(modelFolder)
    except OSError:
        pass

    print(f'Extrapolate {modelName} model results...')

    remapFolder = os.path.join(modelFolder, 'remap')
    try:
        os.makedirs(remapFolder)
    except OSError:
        pass

    remap.remap_model(config, remapFolder)

    _extrap_model(config, modelFolder)

    print('  Done.')


def _extrap_model(config, modelFolder):

    resExtrap = get_res(config, extrap=True)
    resFinal = get_res(config, extrap=False)
    hres = get_horiz_res(config)
    modelName = config.get('model', 'name')

    fields = config.get('model', 'fields')
    fields = fields.replace(',', ' ').split()

    basin = config.get('model', 'basin')
    combineBasins = config.getboolean('model', 'combineBasins')

    inFileName = f'{modelFolder}/remap/{modelName}_temperature_{resExtrap}.nc'
    bedMaskFileName = f'{modelFolder}/bed_mask_{resExtrap}.nc'
    bedFileName = f'bedmap2/bedmap2_{hres}.nc'
    basinNumberFileName = f'imbie/basinNumbers_{hres}.nc'

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    matrixDir = os.path.join(modelName.lower(), 'matrices')
    progressDirs = dict()
    for fieldName in fields:
        progressDirs[fieldName] = f'{modelFolder}/progress_{fieldName}'

    for fieldName in fields:
        progressDir = progressDirs[fieldName]
        inFileName = f'{modelFolder}/remap/' \
                     f'{modelName}_{fieldName}_{resExtrap}.nc'
        outFileName = f'{progressDir}/' \
                      f'{modelName}_{fieldName}_{resExtrap}_extrap_horiz.nc'

        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir, basin=basin, combine=combineBasins)

    if not combineBasins:
        return

    for fieldName in fields:
        progressDir = progressDirs[fieldName]
        inFileName = f'{progressDir}/' \
                     f'{modelName}_{fieldName}_{resExtrap}_extrap_horiz.nc'

        outFileName = f'{progressDir}/' \
                      f'{modelName}_{fieldName}_{resExtrap}_extrap_vert.nc'

        extrap_vert(config, inFileName, outFileName, fieldName)

    inFileNames = {}
    outFileNames = {}
    for fieldName in fields:
        progressDir = progressDirs[fieldName]

        inFileNames[fieldName] = \
            f'{progressDir}/' \
            f'{modelName}_{fieldName}_{resExtrap}_extrap_vert.nc'

        outFileNames[fieldName] = \
            f'{progressDir}/' \
            f'{modelName}_{fieldName}_{resFinal}_extrap_vert.nc'

    remap_vertical(config, inFileNames, outFileNames, extrap=False)

    for fieldName in fields:
        progressDir = progressDirs[fieldName]
        inFileName = f'{progressDir}/' \
                     f'{modelName}_{fieldName}_{resFinal}_extrap_vert.nc'

        outFileName = f'{modelFolder}/{modelName}_{fieldName}_{resFinal}.nc'

        extrap_grounded_above_sea_level(config, inFileName, outFileName,
                                        fieldName, progressDir, matrixDir)
