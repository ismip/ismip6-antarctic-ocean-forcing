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

    inFileName = 'obs/obs_temperature_{}_{}.nc'.format(decades, resExtrap)
    bedMaskFileName = 'obs/bed_mask_{}.nc'.format(resExtrap)
    bedFileName = 'bedmap2/bedmap2_{}.nc'.format(hres)
    basinNumberFileName = 'imbie/basinNumbers_{}.nc'.format(hres)

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    for fieldName in ['temperature', 'salinity']:
        inFileName = 'obs/obs_{}_{}_{}.nc'.format(fieldName, decades,
                                                  resExtrap)
        outFileName = 'obs/obs_{}_{}_{}_extrap_horiz.nc'.format(
                fieldName, decades, resExtrap)

        progressDir = 'obs/progress_{}'.format(fieldName)
        matrixDir = 'obs/matrices_{}'.format(fieldName)
        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir)

    for fieldName in ['temperature', 'salinity']:
        inFileName = 'obs/obs_{}_{}_{}_extrap_horiz.nc'.format(
                fieldName, decades, resExtrap)

        outFileName = 'obs/obs_{}_{}_{}_extrap_vert.nc'.format(
                fieldName, decades, resExtrap)

        extrap_vert(config, inFileName, outFileName, fieldName)

    tempFileName = \
        'obs/obs_temperature_{}_{}_extrap_vert.nc'.format(decades, resExtrap)
    salinFileName = \
        'obs/obs_salinity_{}_{}_extrap_vert.nc'.format(decades, resExtrap)
    outFileName = \
        'obs/obs_thermal_forcing_{}_{}_extrap_vert.nc'.format(decades,
                                                              resExtrap)
    compute_thermal_forcing(tempFileName, salinFileName, outFileName)

    inFileNames = {}
    outFileNames = {}
    for fieldName in ['temperature', 'salinity']:

        inFileNames[fieldName] = \
            'obs/obs_{}_{}_{}_extrap_vert.nc'.format(
                fieldName, decades, resExtrap)

        outFileNames[fieldName] = \
            'obs/obs_{}_{}_{}_extrap_vert.nc'.format(
                fieldName, decades, resFinal)

    remap_vertical(config, inFileNames, outFileNames, extrap=False)

    for fieldName in ['temperature', 'salinity']:

        inFileName = 'obs/obs_{}_{}_{}_extrap_vert.nc'.format(
                fieldName, decades, resFinal)

        outFileName = 'obs/obs_{}_{}_{}.nc'.format(fieldName, decades,
                                                   resFinal)

        progressDir = 'obs/progress_{}'.format(fieldName)
        matrixDir = 'obs/matrices_{}'.format(fieldName)
        extrap_grounded_above_sea_level(config, inFileName, outFileName,
                                        fieldName, progressDir, matrixDir)
