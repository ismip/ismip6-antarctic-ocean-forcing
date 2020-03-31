import xarray
import os
import numpy
import progressbar


def extrap_vert(config, inFileName, outFileName, fieldName, timeIndices=None):

    if os.path.exists(outFileName):
        return

    ds = xarray.open_dataset(inFileName)
    if 'time' in ds.dims and timeIndices is not None:
        ds = ds.isel(time=timeIndices)

    nz = ds.sizes['z']

    field3D = ds[fieldName].values
    origShape = field3D.shape

    if 'time' in ds.dims:
        nt = ds.sizes['time']
    else:
        nt = 1
        field3D = field3D.reshape((1, origShape[0], origShape[1],
                                   origShape[2]))

    print('  Extrapolating {} vertically...'.format(fieldName))
    widgets = ['  z=1/{}: '.format(nz),
               progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nz*nt).start()

    for zIndex in range(1, nz):
        for tIndex in range(nt):

            fieldAbove = field3D[tIndex, zIndex-1, :, :]
            fieldLocal = field3D[tIndex, zIndex, :, :]
            invalid = numpy.isnan(fieldLocal)
            fieldLocal[invalid] = fieldAbove[invalid]
            field3D[tIndex, zIndex, :, :] = fieldLocal

            bar.widgets[0] = '  z={}/{}: '.format(zIndex+1, nz)
            bar.update(tIndex + nt*zIndex)

    bar.finish()

    if 'time' not in ds.dims:
        field3D = field3D.reshape(origShape)

    dims = ds[fieldName].dims
    attrs = ds[fieldName].attrs
    ds[fieldName] = (dims, field3D)
    ds[fieldName].attrs = attrs

    ds.to_netcdf(outFileName)
