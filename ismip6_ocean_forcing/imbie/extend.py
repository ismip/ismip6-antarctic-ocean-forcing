import numpy
import xarray
import imageio
import skfmm
import os


def extend_imbie_masks(res, basins, bedFileName):

    outFileName = 'imbie/basinNumbers_{}.nc'.format(res)
    if os.path.exists(outFileName):
        return

    ds = xarray.Dataset()
    print('  Extending IMBIE basins beyond ice fronts...')
    with xarray.open_dataset(bedFileName) as dsIn:
        nx = dsIn.sizes['x']
        ny = dsIn.sizes['y']
        ds['x'] = dsIn.x
        ds['y'] = dsIn.y

    minDistance = 1e30*numpy.ones((ny, nx))
    basinNumber = -1*numpy.zeros((ny, nx), int)
    for index, basinName in enumerate(basins.keys()):
        print('    {}'.format(basinName))
        imageFileName = 'imbie/basins_{}/{}.png'.format(res, basinName)
        image = imageio.imread(imageFileName)
        basinFraction = 1. - image[::-1, :, 0]/255.
        distance = skfmm.distance(-2.*basinFraction + 1)
        mask = distance < minDistance
        basinNumber[mask] = index
        minDistance[mask] = distance[mask]

    ds['basinNumber'] = (('y', 'x'), basinNumber)
    ds.to_netcdf(outFileName)
