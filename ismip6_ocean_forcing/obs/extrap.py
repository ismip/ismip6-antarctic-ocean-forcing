from ismip6_ocean_forcing.extrap.horiz import make_3D_bed_mask, extrap_horiz, \
    extrap_grounded_above_sea_level
from ismip6_ocean_forcing.extrap.vert import extrap_vert
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res
from ismip6_ocean_forcing.remap.interp1d import remap_vertical
from ismip6_ocean_forcing.thermal_forcing.main import compute_thermal_forcing


def extrap_obs(config, decades):

    resExtrap = get_res(config, extrap=True)
    resFinal = get_res(config, extrap=False)
    hres = get_horiz_res(config)

    inFileName = f'obs/obs_temperature_{decades}_{resExtrap}.nc'
    bedMaskFileName = f'obs/bed_mask_{resExtrap}.nc'
    bedFileName = f'bedmap2/bedmap2_{hres}.nc'
    basinNumberFileName = f'imbie/basinNumbers_{hres}.nc'

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    matrixDirs = dict()
    progressDirs = dict()
    for fieldName in ['temperature', 'salinity']:
        progressDirs[fieldName] = f'obs/progress_{fieldName}'
        matrixDirs[fieldName] = f'obs/matrices_{fieldName}'

    for fieldName in ['temperature', 'salinity']:
        progressDir = progressDirs[fieldName]
        matrixDir = matrixDirs[fieldName]
        inFileName = f'{progressDir}/obs_{fieldName}_{decades}_{resExtrap}.nc'
        outFileName = f'{progressDir}/' \
                      f'obs_{fieldName}_{decades}_{resExtrap}_extrap_horiz.nc'

        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir)

    for fieldName in ['temperature', 'salinity']:
        progressDir = progressDirs[fieldName]
        inFileName = f'{progressDir}/' \
                     f'obs_{fieldName}_{decades}_{resExtrap}_extrap_horiz.nc'
        outFileName = f'{progressDir}/' \
                      f'obs_{fieldName}_{decades}_{resExtrap}_extrap_vert.nc'

        extrap_vert(config, inFileName, outFileName, fieldName)

    tempFileName = \
        f'{progressDir}/obs_temperature_{decades}_{resExtrap}_extrap_vert.nc'
    salinFileName = \
        f'{progressDir}/obs_salinity_{decades}_{resExtrap}_extrap_vert.nc'
    outFileName = \
        f'{progressDir}/' \
        f'obs_thermal_forcing_{decades}_{resExtrap}_extrap_vert.nc'
    compute_thermal_forcing(tempFileName, salinFileName, outFileName)

    inFileNames = {}
    outFileNames = {}
    for fieldName in ['temperature', 'salinity']:
        progressDir = progressDirs[fieldName]

        inFileNames[fieldName] = \
            f'{progressDir}/' \
            f'obs_{fieldName}_{decades}_{resExtrap}_extrap_vert.nc'

        outFileNames[fieldName] = \
            f'{progressDir}/' \
            f'obs_{fieldName}_{decades}_{resFinal}_extrap_vert.nc'

    remap_vertical(config, inFileNames, outFileNames, extrap=False)

    for fieldName in ['temperature', 'salinity']:
        progressDir = progressDirs[fieldName]
        matrixDir = matrixDirs[fieldName]

        inFileName = f'{progressDir}/' \
                     f'obs_{fieldName}_{decades}_{resFinal}_extrap_vert.nc'

        outFileName = f'obs/obs_{fieldName}_{decades}_{resFinal}.nc'

        extrap_grounded_above_sea_level(config, inFileName, outFileName,
                                        fieldName, progressDir, matrixDir)
