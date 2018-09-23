import xarray
import numpy
import os
import zipfile

from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.remap.interp1d import weights_and_indices, interp2d
from ismip6_ocean_forcing.remap.res import get_res


def bedmap2_to_ismip6_grid(config):
    try:
        os.makedirs('bedmap2')
    except OSError:
        pass
    try:
        os.makedirs('ismip6')
    except OSError:
        pass

    res = get_res(config)

    inFileName = 'bedmap2/bedmap2.nc'
    outGridFileName = 'ismip6/{}_grid.nc'.format(res)
    outFileName = 'bedmap2/bedmap2_{}.nc'.format(res)

    _bedmap2_bin_to_netcdf(inFileName)

    _write_ismip6_grid(config, outGridFileName)

    _remap(inFileName, outGridFileName, outFileName, res)


def _write_ismip6_grid(config, outFileName):
    if os.path.exists(outFileName):
        return

    ds = xarray.Dataset()
    dx = config.getfloat('grid', 'dx')
    nx = config.getint('grid', 'nx')
    dy = config.getfloat('grid', 'dy')
    ny = config.getint('grid', 'ny')
    x = dx*numpy.arange(-(nx-1)//2, (nx-1)//2 + 1)
    y = dy*numpy.arange(-(ny-1)//2, (ny-1)//2 + 1)
    ds['x'] = ('x', x)
    ds.x.attrs['units'] = 'meters'
    ds['y'] = ('y', y)
    ds.y.attrs['units'] = 'meters'
    ds.attrs['Grid'] = "Datum = WGS84, earth_radius = 6378137., " \
                       "earth_eccentricity = 0.081819190842621, " \
                       "falseeasting = -3040000., " \
                       "falsenorthing = -3040000., " \
                       "standard_parallel = -71., central_meridien = 0, " \
                       "EPSG=3031"
    ds.attrs['proj'] = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 " \
                       "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    ds.attrs['proj4'] = "+init=epsg:3031"
    ds.to_netcdf(outFileName)


def _bedmap2_bin_to_netcdf(outFileName):

    if os.path.exists(outFileName):
        return

    fields = ['bed', 'surface', 'thickness', 'coverage', 'rockmask',
              'grounded_bed_uncertainty', 'icemask_grounded_and_shelves']

    allExist = True
    for field in fields:
        fileName = 'bedmap2/bedmap2_bin/bedmap2_{}.flt'.format(field)
        if not os.path.exists(fileName):
            allExist = False
            break

    if not allExist:
        # download
        baseURL = 'https://secure.antarctica.ac.uk/data/bedmap2'
        fileNames = ['bedmap2_bin.zip']

        download_files(fileNames, baseURL, 'bedmap2')

        print('Decompressing Bedmap2 data...')
        # unzip
        with zipfile.ZipFile('bedmap2/bedmap2_bin.zip', 'r') as f:
            f.extractall('bedmap2/')
        print('  Done.')

    print('Converting Bedmap2 to NetCDF...')
    ds = xarray.Dataset()
    x = numpy.linspace(-3333000., 3333000., 6667)
    y = x
    ds['x'] = ('x', x)
    ds.x.attrs['units'] = 'meters'
    ds['y'] = ('y', y)
    ds.y.attrs['units'] = 'meters'
    ds.attrs['Grid'] = "Datum = WGS84, earth_radius = 6378137., " \
                       "earth_eccentricity = 0.081819190842621, " \
                       "falseeasting = -3333000., " \
                       "falsenorthing = -3333000., " \
                       "standard_parallel = -71., central_meridien = 0, " \
                       "EPSG=3031"
    ds.attrs['proj'] = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 " \
                       "+x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    ds.attrs['proj4'] = "+init=epsg:3031"

    # add Bedmap2 data
    for fieldName in fields:
        fileName = 'bedmap2/bedmap2_bin/bedmap2_{}.flt'.format(fieldName)
        with open(fileName, 'r') as f:
            field = numpy.fromfile(f, dtype=numpy.float32).reshape(6667, 6667)
            # flip the y axis
            field = field[::-1, :]
            # switch invalid values to be NaN (as expected by xarray)
            field[field == -9999.] = numpy.nan
        if fieldName == 'rockmask':
            # rock mask is zero where rock and -9999 (now NaN) elsewhere
            field = numpy.array(numpy.isfinite(field), numpy.float32)
        if fieldName == 'icemask_grounded_and_shelves':
            # split into separate grounded and floating masks
            ds['icemask_grounded'] = \
                (('y', 'x'), numpy.array(field == 0, numpy.float32))
            ds['icemask_shelves'] = \
                (('y', 'x'), numpy.array(field == 1, numpy.float32))
            ds['open_ocean_mask'] = \
                (('y', 'x'), numpy.array(numpy.isnan(field), numpy.float32))
        else:
            ds[fieldName] = (('y', 'x'), field)

    ds.to_netcdf(outFileName)
    print('  Done.')


def _remap(inFileName, outGridFileName, outFileName, res):

    if os.path.exists(outFileName):
        return

    print('Remapping Bedmap2 to {} grid...'.format(res))
    dsIn = xarray.open_dataset(inFileName)
    dsOut = xarray.open_dataset(outGridFileName)

    print('  Computing remapping weights...')
    xWeights, xIndices = weights_and_indices(xInCenter=dsIn.x.values,
                                             xOutCenter=dsOut.x.values,
                                             xDim='xOut')

    yWeights, yIndices = weights_and_indices(xInCenter=dsIn.y.values,
                                             xOutCenter=dsOut.y.values,
                                             xDim='yOut')

    fields = ['bed', 'surface', 'thickness', 'rockmask',
              'grounded_bed_uncertainty', 'icemask_grounded',
              'icemask_shelves', 'open_ocean_mask']

    print('  Remapping fields...')
    for fieldName in fields:
        print('    {}'.format(fieldName))
        result = interp2d(dsIn[fieldName], xWeights, xIndices,
                          yWeights, yIndices, normalizationThreshold=0.1)
        result = result.rename({'xOut': 'x', 'yOut': 'y'})
        dsOut[fieldName] = result
    dsOut.to_netcdf(outFileName)
    print('  Done.')
