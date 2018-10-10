from ismip6_ocean_forcing.extrap.horiz import make_3D_bed_mask, extrap_horiz, \
    extrap_grounded_above_sea_level
from ismip6_ocean_forcing.extrap.vert import extrap_vert
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res
from ismip6_ocean_forcing.remap.interp1d import remap_vertical
from ismip6_ocean_forcing.thermal_forcing.main import compute_thermal_forcing


def extrap_woa(config):

    resExtrap = get_res(config, extrap=True)
    resFinal = get_res(config, extrap=False)
    hres = get_horiz_res(config)

    inFileName = 'woa/woa_temperature_1955-2012_{}.nc'.format(resExtrap)
    bedMaskFileName = 'woa/bed_mask_{}.nc'.format(resExtrap)
    bedFileName = 'bedmap2/bedmap2_{}.nc'.format(hres)
    basinNumberFileName = 'imbie/basinNumbers_{}.nc'.format(hres)

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    for fieldName in ['temperature', 'salinity']:
        inFileName = 'woa/woa_{}_1955-2012_{}.nc'.format(fieldName, resExtrap)
        outFileName = 'woa/woa_{}_1955-2012_{}_extrap_horiz.nc'.format(
                fieldName, resExtrap)

        progressDir = 'woa/progress_{}'.format(fieldName)
        matrixDir = 'woa/matrices_{}'.format(fieldName)
        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir)

    for fieldName in ['temperature', 'salinity']:
        inFileName = 'woa/woa_{}_1955-2012_{}_extrap_horiz.nc'.format(
                fieldName, resExtrap)

        outFileName = 'woa/woa_{}_1955-2012_{}_extrap_vert.nc'.format(
                fieldName, resExtrap)

        extrap_vert(config, inFileName, outFileName, fieldName)

    tempFileName = \
        'woa/woa_temperature_1955-2012_{}_extrap_vert.nc'.format(resExtrap)
    salinFileName = \
        'woa/woa_salinity_1955-2012_{}_extrap_vert.nc'.format(resExtrap)
    outFileName = \
        'woa/woa_thermal_forcing_1955-2012_{}_extrap_vert.nc'.format(resExtrap)
    compute_thermal_forcing(tempFileName, salinFileName, outFileName)

    inFileNames = {}
    outFileNames = {}
    for fieldName in ['temperature', 'salinity']:

        inFileNames[fieldName] = \
            'woa/woa_{}_1955-2012_{}_extrap_vert.nc'.format(
                fieldName, resExtrap)

        outFileNames[fieldName] = \
            'woa/woa_{}_1955-2012_{}_extrap_vert.nc'.format(
                fieldName, resFinal)

    remap_vertical(config, inFileNames, outFileNames, extrap=False)

    for fieldName in ['temperature', 'salinity']:

        inFileName = 'woa/woa_{}_1955-2012_{}_extrap_vert.nc'.format(
                fieldName, resFinal)

        outFileName = 'woa/woa_{}_1955-2012_{}.nc'.format(fieldName, resFinal)

        progressDir = 'woa/progress_{}'.format(fieldName)
        matrixDir = 'woa/matrices_{}'.format(fieldName)
        extrap_grounded_above_sea_level(config, inFileName, outFileName,
                                        fieldName, progressDir, matrixDir)
