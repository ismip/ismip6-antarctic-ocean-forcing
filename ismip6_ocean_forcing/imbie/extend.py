import numpy
import xarray
from scipy import misc
import skfmm
import os


def extend_imbie_masks(basins, bedFileName):

    outFileName = 'imbie/basinNumbers.nc'
    if os.path.exists(outFileName):
        return

    ds = xarray.Dataset()
    print('  Extending IMBIE basins beyond ice fronts...')
    with xarray.open_dataset(bedFileName) as dsIn:
        nx = dsIn.sizes['x']
        ny = dsIn.sizes['y']
        ds['x'] = dsIn.x
        ds['y'] = dsIn.y

    basinNames = []
    for basinList in basins:
        if isinstance(basinList, list):
            basinName = 'Antarctica_IMBIE{}'.format(
                    '_'.join(['{}'.format(basin) for basin in basinList]))
            basinNames.append(basinName)

        else:
            basinNames.append('Antarctica_IMBIE{}'.format(basinList))

    minDistance = 1e30*numpy.ones((ny, nx))
    basinNumber = -1*numpy.zeros((ny, nx), int)
    for index, basinName in enumerate(basinNames):
        print('    {}'.format(basinName))
        imageFileName = 'imbie/basins/{}.png'.format(basinName)
        image = misc.imread(imageFileName)
        basinFraction = 1. - image[::-1, :, 0]/255.
        distance = skfmm.distance(-2.*basinFraction + 1)
        mask = distance < minDistance
        basinNumber[mask] = index
        minDistance[mask] = distance[mask]

    ds['basinNumber'] = (('y', 'x'), basinNumber)
    ds.to_netcdf(outFileName)
