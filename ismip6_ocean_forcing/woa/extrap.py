from ismip6_ocean_forcing.extrap.horiz import make_3D_bed_mask, extrap_horiz
from ismip6_ocean_forcing.extrap.vert import extrap_vert
from ismip6_ocean_forcing.remap.res import get_res


def extrap_woa(config):

    res = get_res(config)
    inFileName = 'woa/woa_temperature_1995-2012_{}.nc'.format(res)
    bedMaskFileName = 'woa/bed_mask_{}.nc'.format(res)
    bedFileName = 'bedmap2/bedmap2_{}.nc'.format(res)
    basinNumberFileName = 'imbie/basinNumbers.nc'

    make_3D_bed_mask(inFileName, bedMaskFileName, bedFileName)

    for fieldName in ['temperature', 'salinity']:
        inFileName = 'woa/woa_{}_1995-2012_{}.nc'.format(fieldName, res)
        outFileName = \
            'woa/woa_{}_1995-2012_{}_extrap_horiz.nc'.format(fieldName, res)

        progressDir = 'woa/progress_{}'.format(fieldName)
        matrixDir = 'woa/matrices'
        extrap_horiz(config, inFileName, outFileName, fieldName, bedFileName,
                     basinNumberFileName, bedMaskFileName, progressDir,
                     matrixDir)

    for fieldName in ['temperature', 'salinity']:
        inFileName = \
            'woa/woa_{}_1995-2012_{}_extrap_horiz.nc'.format(fieldName, res)

        outFileName = \
            'woa/woa_{}_1995-2012_{}_extrap_vert.nc'.format(fieldName, res)

        extrap_vert(config, inFileName, outFileName, fieldName)
