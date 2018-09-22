import xarray
import numpy
import os
from scipy.signal import convolve2d
import skfmm
import progressbar


def make_3D_bed_mask(inFileName, outFileName, bedFileName):

    if os.path.exists(outFileName):
        return

    with xarray.open_dataset(bedFileName) as ds:
        bed = ds.bed.values

    dsOut = xarray.Dataset()
    with xarray.open_dataset(inFileName) as ds:
        nx = ds.sizes['x']
        ny = ds.sizes['y']
        nz = ds.sizes['z']
        for var in ['x', 'y', 'z', 'z_bnds']:
            dsOut[var] = ds[var]

        zTop = ds.z_bnds[:, 0].values

        bedMask = numpy.zeros((nz, nx, ny), bool)

        for zIndex in range(nz):
            mask = numpy.isfinite(bed)
            mask[mask] = bed[mask] <= zTop[zIndex]
            bedMask[zIndex, :, :] = mask

        dsOut['bedMask'] = (('z', 'x', 'y'), bedMask)
        dsOut.to_netcdf(outFileName)


def extrap_horiz(config, inFileName, outFileName, fieldName, bedmap2FileName,
                 basinNumberFileName, bedMaskFileName, progressDir):

    if os.path.exists(outFileName):
        return

    with xarray.open_dataset(bedmap2FileName) as dsBed:
        openOceanMask = dsBed.open_ocean_mask.values > 0.

    with xarray.open_dataset(basinNumberFileName) as dsBasin:
        basinNumbers = dsBasin.basinNumber.values

    try:
        os.makedirs(progressDir)
    except OSError:
        pass

    progressFileName = '{}/open_ocean.nc'.format(progressDir)
    if os.path.exists(progressFileName):
        ds = xarray.open_dataset(progressFileName)
    else:
        ds = xarray.open_dataset(inFileName)

        # mask out bed and areas under ice shelves or in grounded ice regions
        _mask_ice_and_bed(ds, fieldName, openOceanMask, bedMaskFileName)

        # first, extrapolate the open ocean
        validMask = openOceanMask
        invalidMask = openOceanMask
        _extrap_basin(config, ds, fieldName, 'open ocean', validMask,
                      invalidMask,  bedMaskFileName)

        ds.to_netcdf(progressFileName)

    # then, extrapolate each IMBIE basin
    basinCount = numpy.amax(basinNumbers) + 1
    for basinNumber in range(basinCount):

        progressFileName = '{}/basin{}.nc'.format(progressDir, basinNumber)
        if os.path.exists(progressFileName):
            ds = xarray.open_dataset(progressFileName)
        else:
            validMask = openOceanMask
            invalidMask = numpy.logical_or(basinNumbers == basinNumber,
                                           openOceanMask)
            basinName = 'basin {}/{}'.format(basinNumber+1, basinCount)
            _extrap_basin(config, ds, fieldName, basinName, validMask,
                          invalidMask, bedMaskFileName)
            ds.to_netcdf(progressFileName)

    ds.to_netcdf(outFileName)


def _mask_ice_and_bed(ds, fieldName, openOceanMask, bedMaskFileName):

    dsMask = xarray.open_dataset(bedMaskFileName)

    field3D = ds[fieldName].values
    for zIndex in range(ds.sizes['z']):
        bedMask = dsMask.bedMask[zIndex, :, :].values
        mask = numpy.logical_not(numpy.logical_and(bedMask, openOceanMask))
        if 'time' in ds.dims:
            for tIndex in range(ds.sizes['time']):
                field = field3D[tIndex, zIndex, :, :]
                field[mask] = numpy.NaN
                field3D[tIndex, zIndex, :, :] = field
        else:
            field = field3D[zIndex, :, :]
            field[mask] = numpy.NaN
            field3D[zIndex, :, :] = field

    dims = ds[fieldName].dims
    attrs = ds[fieldName].attrs
    ds[fieldName] = (dims, field3D)
    ds[fieldName].attrs = attrs


def _extrap_basin(config, ds, fieldName, basinName, validMask, invalidMask,
                  bedMaskFileName):

    dsMask = xarray.open_dataset(bedMaskFileName)

    nz = ds.sizes['z']

    dx = config.getfloat('grid', 'dx')
    kernelRadius = config.getfloat('extrapolation', 'kernelRadius')
    maxDistance = config.getfloat('extrapolation', 'maxDistance')
    weightSumThreshold = config.getfloat('extrapolation', 'weightSumThreshold')

    # the kernel should be big enough to capture weights up to 0.01 of the peak
    kernelSize = int(numpy.ceil(kernelRadius*3/dx))
    x = dx*numpy.arange(-kernelSize, kernelSize+1)/kernelRadius
    x, y = numpy.meshgrid(x, x)

    kernel = numpy.exp(-0.5*(x**2 + y**2))
    kernel = kernel/numpy.sum(kernel)

    field3D = ds[fieldName].values

    origShape = field3D.shape

    if 'time' in ds.dims:
        nt = ds.sizes['time']
    else:
        nt = 1
        field3D = field3D.reshape((1, origShape[0], origShape[1],
                                   origShape[2]))

    workCount = 0
    for zIndex in range(nz):
        bedMask = dsMask.bedMask[zIndex, :, :].values
        field = field3D[0, zIndex, :, :]
        valid = numpy.logical_and(numpy.isfinite(field), validMask)
        invalid = numpy.logical_and(numpy.isnan(field), invalidMask)
        fillMask = numpy.logical_and(invalid, bedMask)
        invalidCount = numpy.count_nonzero(fillMask)
        workCount += invalidCount

    print('  Extrapolating {} in {}...'.format(fieldName, basinName))
    widgets = ['  z=1/{}: '.format(nz),
               progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=workCount).start()

    progress = 0
    for zIndex in range(nz):

        bedMask = dsMask.bedMask[zIndex, :, :].values
        field = field3D[0, zIndex, :, :]
        valid = numpy.logical_and(numpy.isfinite(field), validMask)
        invalid = numpy.logical_and(numpy.isnan(field), invalidMask)
        fillMask = numpy.logical_and(invalid, bedMask)
        invalidCount = numpy.count_nonzero(fillMask)
        localWork = invalidCount

        while True:

            distance = skfmm.distance(-2.*valid + 1., dx=dx)

            maskExtrap = convolve2d(valid, kernel, mode='same')
            mask = numpy.logical_and(fillMask, maskExtrap > weightSumThreshold)
            mask = numpy.logical_and(mask, distance <= maxDistance)

            newValid = numpy.logical_and(numpy.logical_not(valid), mask)
            newValidCount = numpy.count_nonzero(newValid)
            invalidCount -= newValidCount
            # print(zIndex, newValidCount, invalidCount)
            if newValidCount == 0:
                break

            for tIndex in range(nt):
                fieldExtrap = field3D[tIndex, zIndex, :, :].copy()

                fieldExtrap[numpy.logical_not(valid)] = 0.
                fieldExtrap[numpy.isnan(fieldExtrap)] = 0.
                fieldExtrap = convolve2d(fieldExtrap, kernel, mode='same')

                field[mask] = fieldExtrap[mask]/maskExtrap[mask]
                field3D[tIndex, zIndex, :, :] = field

            valid[mask] = True
            widgets[0] = '  z={}/{}: '.format(zIndex+1, nz)
            bar.update(progress + localWork - invalidCount)

        progress += localWork

    bar.finish()

    if 'time' not in ds.dims:
        field3D = field3D.reshape(origShape)

    dims = ds[fieldName].dims
    attrs = ds[fieldName].attrs
    ds[fieldName] = (dims, field3D)
    ds[fieldName].attrs = attrs
