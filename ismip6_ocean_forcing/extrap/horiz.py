import xarray
import numpy
import os
from scipy.signal import convolve2d
from scipy.ndimage.morphology import binary_fill_holes
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
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
                 basinNumberFileName, bedMaskFileName, progressDir,
                 matrixDir):

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

    try:
        os.makedirs(matrixDir)
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
        basinName = 'open ocean'

        matrixFileTemplate = '{}/matrix_open_ocean_{{}}.npz'.format(matrixDir)
        _write_basin_matrices(config, ds, fieldName, basinName, validMask,
                              invalidMask,  bedMaskFileName,
                              matrixFileTemplate)
        _extrap_basin(config, ds, fieldName, basinName, matrixFileTemplate)

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

            matrixFileTemplate = '{}/matrix_basin{}_{{}}.npz'.format(
                    matrixDir, basinNumber)
            _write_basin_matrices(config, ds, fieldName, basinName, validMask,
                                  invalidMask,  bedMaskFileName,
                                  matrixFileTemplate)

            _extrap_basin(config, ds, fieldName, basinName, matrixFileTemplate)
            ds.to_netcdf(progressFileName)

    ds.to_netcdf(outFileName)


def _mask_ice_and_bed(ds, fieldName, openOceanMask, bedMaskFileName):

    openOceanMask = numpy.logical_and(openOceanMask, ds.lat.values < -60.)

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


def _write_basin_matrices(config, ds, fieldName, basinName, validMask,
                          invalidMask, bedMaskFileName, matrixFileTemplate):

    dsMask = xarray.open_dataset(bedMaskFileName)

    validMask = numpy.logical_and(validMask, ds.lat.values < -60.)
    invalidMask = numpy.logical_and(invalidMask, ds.lat.values < -60.)

    nz = ds.sizes['z']

    dx = config.getfloat('grid', 'dx')
    kernelRadius = config.getfloat('extrapolation', 'kernelRadius')

    # the kernel should be big enough to capture weights up to 0.01 of the peak
    kernelSize = int(numpy.ceil(kernelRadius*3/dx))
    x = dx*numpy.arange(-kernelSize, kernelSize+1)/kernelRadius
    x, y = numpy.meshgrid(x, x)

    kernel = numpy.exp(-0.5*(x**2 + y**2))
    kernel = kernel/numpy.sum(kernel)

    field3D = ds[fieldName].values

    if 'time' in ds.dims:
        field3D = field3D[0, :, :, :]

    print('  Writing matrices for {} in {}...'.format(fieldName, basinName))

    widgets = ['  ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nz).start()

    for zIndex in range(nz):

        outFileName = matrixFileTemplate.format(zIndex)
        if os.path.exists(outFileName):
            continue

        bedMask = dsMask.bedMask[zIndex, :, :].values
        field = field3D[zIndex, :, :]
        valid = numpy.logical_and(numpy.isfinite(field), validMask)
        fillMask = numpy.logical_and(numpy.isnan(field), invalidMask)
        fillMask = numpy.logical_and(fillMask, bedMask)
        dataMask = numpy.logical_or(valid, fillMask)
        dataMask = numpy.logical_not(binary_fill_holes(
                numpy.logical_not(dataMask)))
        valid = numpy.logical_and(valid, dataMask)
        fillMask = numpy.logical_and(fillMask, dataMask)
        fillCount = numpy.count_nonzero(fillMask)

        weightSum = convolve2d(dataMask, kernel, mode='same')

        ny, nx = fillMask.shape
        indices = numpy.indices((ny, nx))
        yIndices = indices[0].ravel()
        xIndices = indices[1].ravel()

        fillIndices = numpy.nonzero(fillMask.ravel())[0]
        weightSum = weightSum[fillMask]

        fillInverseMap = -1*numpy.ones((ny, nx), int)
        fillInverseMap[fillMask] = numpy.arange(fillCount)

        matrix = lil_matrix((fillCount, fillCount))
        for index, fillIndex in enumerate(fillIndices):
            xc = xIndices[fillIndex]
            yc = yIndices[fillIndex]
            xMin = max(0, xc-kernelSize)
            xMax = min(nx, xc+kernelSize+1)
            yMin = max(0, yc-kernelSize)
            yMax = min(ny, yc+kernelSize+1)
            kxMin = xMin-xc+kernelSize
            kxMax = xMax-xc+kernelSize
            kyMin = yMin-yc+kernelSize
            kyMax = yMax-yc+kernelSize
            otherIndices = fillInverseMap[yMin:yMax, xMin:xMax]
            weights = kernel[kyMin:kyMax, kxMin:kxMax]/weightSum[index]
            mask = dataMask[yMin:yMax, xMin:xMax]
            mask = otherIndices >= 0
            otherIndices = otherIndices[mask]
            weights = weights[mask]
            matrix[index, otherIndices] = -weights
            # put ones along the diagonal
            matrix[index, index] = 1 + matrix[index, index]

        matrix = matrix.tocsr()
        _save_matrix_and_kernel(outFileName, matrix, kernel, valid, fillMask,
                                weightSum)
        bar.update(zIndex)
    bar.finish()


def _extrap_basin(config, ds, fieldName, basinName, matrixFileTemplate):

    nz = ds.sizes['z']

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

        matrix, kernel, valid, fillMask, weightSum = _load_matrix_and_kernel(
                matrixFileTemplate.format(zIndex))

        fillCount = numpy.count_nonzero(fillMask)

        workCount += nt*fillCount

    print('  Extrapolating {} in {}...'.format(fieldName, basinName))
    widgets = ['  z=1/{}: '.format(nz),
               progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=workCount).start()

    progress = 0
    for zIndex in range(nz):

        matrix, kernel, valid, fillMask, weightSum = _load_matrix_and_kernel(
                matrixFileTemplate.format(zIndex))

        fillCount = numpy.count_nonzero(fillMask)
        if fillCount == 0:
            continue

        widgets[0] = '  z={}/{}: '.format(zIndex+1, nz)

        for tIndex in range(nt):
            fieldSlice = field3D[tIndex, zIndex, :, :]
            fieldExtrap = fieldSlice.copy()
            fieldExtrap[numpy.logical_not(valid)] = 0.
            fieldExtrap[numpy.isnan(fieldExtrap)] = 0.
            fieldExtrap = convolve2d(fieldExtrap, kernel, mode='same')
            rhs = fieldExtrap[fillMask]/weightSum

            fieldFill = spsolve(matrix, rhs)
            fieldSlice[fillMask] = fieldFill

            field3D[tIndex, zIndex, :, :] = fieldSlice

            progress += fillCount
            bar.update(progress)

    bar.finish()

    if 'time' not in ds.dims:
        field3D = field3D.reshape(origShape)

    dims = ds[fieldName].dims
    attrs = ds[fieldName].attrs
    ds[fieldName] = (dims, field3D)
    ds[fieldName].attrs = attrs


def _save_matrix_and_kernel(fileName, matrix, kernel, valid, fillMask,
                            weightSum):
    numpy.savez(fileName, data=matrix.data, indices=matrix.indices,
                indptr=matrix.indptr, shape=matrix.shape, kernel=kernel,
                valid=valid, fillMask=fillMask, weightSum=weightSum)


def _load_matrix_and_kernel(fileName):
    loader = numpy.load(fileName)
    kernel = loader['kernel']
    valid = loader['valid']
    fillMask = loader['fillMask']
    weightSum = loader['weightSum']
    matrix = csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                        shape=loader['shape'])
    return matrix, kernel,  valid, fillMask, weightSum
