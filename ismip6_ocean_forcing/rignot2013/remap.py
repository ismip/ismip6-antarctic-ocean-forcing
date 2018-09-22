import xarray
import numpy
import os

from ismip6_ocean_forcing.remap.interp1d import weights_and_indices, interp2d
from ismip6_ocean_forcing.remap.res import get_res


def rignot_to_ismip6_grid(config):

    if not config.getboolean('rignot', 'remap'):
        return

    try:
        os.makedirs('rignot')
    except OSError:
        pass

    res = get_res(config)

    inFileName = config.get('rignot', 'fileName')

    outGridFileName = 'ismip6/{}_grid.nc'.format(res)
    outFileName = 'rignot/rignot_melt_rates_{}.nc'.format(res)

    _remap(inFileName, outGridFileName, outFileName, res)


def _remap(inFileName, outGridFileName, outFileName, res):

    if os.path.exists(outFileName):
        return

    print('Remapping Rignot et al. (2013) melt rates to {} grid...'.format(
            res))
    dsIn = xarray.open_dataset(inFileName)
    dsIn = dsIn.rename({'nx': 'x', 'xaxis': 'x', 'ny': 'y', 'yaxis': 'y'})
    ny = dsIn.sizes['y']
    yIndices = numpy.arange(ny-1, -1, -1)
    dsIn = dsIn.isel(y=yIndices)
    dsOut = xarray.open_dataset(outGridFileName)

    print('  Computing remapping weights...')
    xWeights, xIndices = weights_and_indices(xInCenter=dsIn.x.values,
                                             xOutCenter=dsOut.x.values,
                                             xDim='xOut')

    yWeights, yIndices = weights_and_indices(xInCenter=dsIn.y.values,
                                             xOutCenter=dsOut.y.values,
                                             xDim='yOut')

    fields = ['melt_actual', 'melt_steadystate', 'lon', 'lat']

    print('  Remapping fields...')
    for fieldName in fields:
        print('    {}'.format(fieldName))
        result = interp2d(dsIn[fieldName], xWeights, xIndices,
                          yWeights, yIndices, normalizationThreshold=None)
        result = result.rename({'xOut': 'x', 'yOut': 'y'})
        dsOut[fieldName] = result
    dsOut.to_netcdf(outFileName)
    print('  Done.')
