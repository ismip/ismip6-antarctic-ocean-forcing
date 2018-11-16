import xarray
import numpy
import os
from scipy.signal import convolve2d
from scipy.ndimage.morphology import binary_fill_holes
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve
import progressbar
import skfmm
from multiprocessing import Pool
from functools import partial


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
        openOceanMask = dsBed.open_ocean_mask.values >= 0.5
        bed = dsBed.bed.values
        mask = numpy.isfinite(bed)

        # continental shelf mask is all places that have ocean depth less than
        # 1500 m and are connected to the edges of the domain (open ocean)
        continentalShelfMask = numpy.zeros(bed.shape)
        continentalShelfMask[mask] = bed[mask] > -1500.
        # flood fill to take out deep areas inside the continent
        continentalShelfMask = binary_fill_holes(continentalShelfMask)

        # if ice is over ocean deeper than 1500 m, put it in the open ocean
        # instead
        openOceanMask = numpy.logical_or(
                openOceanMask, numpy.logical_not(continentalShelfMask))

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

    validKernelRadius = config.getfloat('extrapolation', 'validKernelRadius')
    invalidKernelRadius = config.getfloat('extrapolation',
                                          'invalidKernelRadius')

    validSmoothingIterations = config.getint('extrapolation',
                                             'validSmoothingIterations')

    dx = config.getfloat('grid', 'dx')
    parallelTasks = config.getint('parallel', 'tasks')

    maskedFileName = '{}/masked.nc'.format(progressDir)
    if not os.path.exists(maskedFileName):
        ds = xarray.open_dataset(inFileName)
        # mask out bed and areas under ice shelves or in grounded ice regions
        _mask_ice_and_bed(ds, fieldName, openOceanMask, bedMaskFileName)
        ds.to_netcdf(maskedFileName)

    # write matrices for each basin and vertical level
    progressFileName = '{}/open_ocean.nc'.format(progressDir)
    if not os.path.exists(progressFileName):
        ds = xarray.open_dataset(maskedFileName)

        # first, extrapolate the open ocean
        validMask = openOceanMask
        invalidMask = openOceanMask
        basinMask = openOceanMask
        basinName = 'open ocean'

        matrixFileTemplate = '{}/matrix_open_ocean_{{}}.npz'.format(matrixDir)
        _write_basin_matrices(ds, fieldName, basinName, openOceanMask,
                              validMask, invalidMask, basinMask,
                              validKernelRadius, invalidKernelRadius,
                              validSmoothingIterations, dx,
                              matrixFileTemplate, parallelTasks,
                              bedMaskFileName)

    basinCount = numpy.amax(basinNumbers) + 1

    basinMaskFileName = '{}/basin_masks.nc'.format(matrixDir)
    if os.path.exists(basinMaskFileName):
        dsBasinMasks = xarray.open_dataset(basinMaskFileName)
    else:
        dsBasinMasks = xarray.Dataset()
        for basinNumber in range(basinCount):
            basinMask = _compute_valid_basin_mask(
                    basinNumbers, basinNumber, openOceanMask,
                    continentalShelfMask, validKernelRadius, dx)
            dsBasinMasks['basin{}Mask'.format(basinNumber)] = \
                (('y', 'x'), basinMask)
        dsBasinMasks.to_netcdf(basinMaskFileName)

    for basinNumber in range(basinCount):

        progressFileName = '{}/basin{}.nc'.format(progressDir, basinNumber)
        if not os.path.exists(progressFileName):
            ds = xarray.open_dataset(maskedFileName)

            basinMask = dsBasinMasks['basin{}Mask'.format(basinNumber)].values
            validMask = basinMask.copy()
            invalidMask = basinMask.copy()
            basinName = 'basin {}/{}'.format(basinNumber+1, basinCount)

            matrixFileTemplate = '{}/matrix_basin{}_{{}}.npz'.format(
                    matrixDir, basinNumber)
            _write_basin_matrices(ds, fieldName, basinName, openOceanMask,
                                  validMask, invalidMask, basinMask,
                                  validKernelRadius, invalidKernelRadius,
                                  validSmoothingIterations, dx,
                                  matrixFileTemplate, parallelTasks,
                                  bedMaskFileName)

    # extrapolate each basin and vertical level
    dsOut = xarray.open_dataset(maskedFileName)

    for basinNumber in range(basinCount):

        basinMask = numpy.logical_and(basinNumbers == basinNumber,
                                      numpy.logical_not(openOceanMask))

        progressFileName = '{}/basin{}.nc'.format(progressDir, basinNumber)
        if os.path.exists(progressFileName):
            ds = xarray.open_dataset(progressFileName)
        else:
            ds = xarray.open_dataset(maskedFileName)

            basinName = 'basin {}/{}'.format(basinNumber+1, basinCount)

            matrixFileTemplate = '{}/matrix_basin{}_{{}}.npz'.format(
                    matrixDir, basinNumber)

            _extrap_basin(ds, fieldName, basinName, matrixFileTemplate,
                          parallelTasks, basinMask, validSmoothingIterations,
                          replaceValidWithSmoothed=True)
            ds.to_netcdf(progressFileName)

        _add_basin_field(ds, dsOut, fieldName, basinMask)

    progressFileName = '{}/open_ocean.nc'.format(progressDir)
    if os.path.exists(progressFileName):
        ds = xarray.open_dataset(progressFileName)
    else:
        ds = xarray.open_dataset(maskedFileName)

        # first, extrapolate the open ocean
        basinName = 'open ocean'

        matrixFileTemplate = '{}/matrix_open_ocean_{{}}.npz'.format(matrixDir)
        _extrap_basin(ds, fieldName, basinName, matrixFileTemplate,
                      parallelTasks, openOceanMask, validSmoothingIterations,
                      replaceValidWithSmoothed=True)

        ds.to_netcdf(progressFileName)

    _add_basin_field(ds, dsOut, fieldName, openOceanMask)

    dsOut.to_netcdf(outFileName)


def extrap_grounded_above_sea_level(config, inFileName, outFileName, fieldName,
                                    progressDir, matrixDir):

    if os.path.exists(outFileName):
        return

    try:
        os.makedirs(progressDir)
    except OSError:
        pass

    try:
        os.makedirs(matrixDir)
    except OSError:
        pass

    ds = xarray.open_dataset(inFileName)

    # first, extrapolate the open ocean
    nx = ds.sizes['x']
    ny = ds.sizes['y']
    openOceanMask = numpy.ones((ny, nx), bool)
    validMask = numpy.ones((ny, nx), bool)
    invalidMask = numpy.ones((ny, nx), bool)
    basinMask = numpy.ones((ny, nx), bool)

    invalidKernelRadius = config.getfloat('extrapolation',
                                          'invalidKernelRadius')
    # we're extrapolating from data that's already been extrapolated, so no
    # no need to use the larger valid radius
    validKernelRadius = invalidKernelRadius
    # ... and no need to do multiple iterations
    validSmoothingIterations = 1

    dx = config.getfloat('grid', 'dx')
    parallelTasks = config.getint('parallel', 'tasks')

    basinName = 'grounded above sea level'

    matrixFileTemplate = '{}/matrix_grounded_above_sea_level_{{}}.npz'.format(
            matrixDir)
    _write_basin_matrices(ds, fieldName, basinName, openOceanMask, validMask,
                          invalidMask, basinMask, validKernelRadius,
                          invalidKernelRadius, validSmoothingIterations,
                          dx, matrixFileTemplate, parallelTasks)

    _extrap_basin(ds, fieldName, basinName, matrixFileTemplate,
                  parallelTasks, basinMask, validSmoothingIterations,
                  replaceValidWithSmoothed=False)

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


def _compute_valid_basin_mask(basinNumbers, basin, openOceanMask,
                              continentalShelfMask, validKernelRadius, dx):

    basinMask = numpy.logical_and(basinNumbers == basin,
                                  numpy.logical_not(openOceanMask))
    otherBasins = numpy.logical_and(basinNumbers != basin,
                                    numpy.logical_not(openOceanMask))
    mask = numpy.logical_or(otherBasins,
                            numpy.logical_not(continentalShelfMask))
    phi = numpy.ma.masked_array(-2.*basinMask + 1., mask=mask)

    distance = skfmm.distance(phi, dx=dx)
    distance = distance.filled(fill_value=0.)
    validMask = numpy.logical_and(distance > 0.,
                                  distance < 3*validKernelRadius)
    basinMask = numpy.logical_or(basinMask, validMask)
    basinMask = numpy.logical_and(basinMask, continentalShelfMask)
    return basinMask


def _write_basin_matrices(ds, fieldName, basinName, openOceanMask, validMask,
                          invalidMask, basinMask, validKernelRadius,
                          invalidKernelRadius, validSmoothingIterations,
                          dx, matrixFileTemplate, parallelTasks,
                          bedMaskFileName=None):

    def get_kernel(kernelRadius):
        # the kernel should be big enough to capture weights up to 0.01 of the
        # peak
        kernelSize = int(numpy.ceil(kernelRadius*3/dx))
        x = dx*numpy.arange(-kernelSize, kernelSize+1)/kernelRadius
        x, y = numpy.meshgrid(x, x)

        kernel = numpy.exp(-0.5*(x**2 + y**2))
        # kernel = kernel/numpy.sum(kernel)
        return kernelSize, kernel

    nz = ds.sizes['z']

    allExist = True
    for zIndex in range(nz):
        if not os.path.exists(matrixFileTemplate.format(zIndex)):
            allExist = False

    if allExist:
        return

    field3D = ds[fieldName].values
    if 'time' in ds.dims:
        field3D = field3D[0, :, :, :]

    if bedMaskFileName is None:
        bedMask = numpy.ones(field3D.shape)
    else:
        with xarray.open_dataset(bedMaskFileName) as dsMask:
            bedMask = dsMask.bedMask.values

    validMask = numpy.logical_and(validMask, ds.lat.values < -60.)
    invalidMask = numpy.logical_and(invalidMask, ds.lat.values < -60.)

    validKernelSize, validKernel = get_kernel(validKernelRadius)
    invalidKernelSize, invalidKernel = get_kernel(invalidKernelRadius)

    print('  Writing matrices for {} in {}...'.format(fieldName, basinName))

    pool = Pool(parallelTasks)
    zIndices = range(nz)
    widgets = ['  ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nz).start()

    partial_func = partial(_write_level_basin_matrix, matrixFileTemplate,
                           field3D, bedMask, validMask, invalidMask,
                           basinMask, openOceanMask, validKernel,
                           invalidKernel, invalidKernelSize,
                           validSmoothingIterations)
    for zIndex, _ in enumerate(pool.imap(partial_func, zIndices)):
        bar.update(zIndex+1)
    bar.finish()


def _write_level_basin_matrix(matrixFileTemplate, field3D, bedMask, validMask,
                              invalidMask, basinMask, openOceanMask,
                              validKernel, invalidKernel, invalidKernelSize,
                              validSmoothingIterations, zIndex):
    outFileName = matrixFileTemplate.format(zIndex)
    if os.path.exists(outFileName):
        return

    field = field3D[zIndex, :, :]

    valid = numpy.logical_and(numpy.isfinite(field), validMask)
    fillMask = numpy.logical_and(numpy.isnan(field), invalidMask)

    dataMask = numpy.logical_or(openOceanMask, fillMask)
    dataMask = numpy.logical_or(dataMask, valid)
    dataMask = numpy.logical_and(dataMask, bedMask[zIndex, :, :])

    # flood fill to make sure the data is contiguous
    dataMask = numpy.logical_not(binary_fill_holes(
            numpy.logical_not(dataMask)))

    # mask to the basin, once we've done the flood fill
    dataMask = numpy.logical_and(dataMask, basinMask)

    # only take valid data and fill data that's contiguous and in the basin
    valid = numpy.logical_and(valid, dataMask)
    fillMask = numpy.logical_and(fillMask, dataMask)
    fillCount = numpy.count_nonzero(fillMask)

    validWeightSum = valid.copy()
    iterMask = numpy.logical_not(dataMask)
    for iterIndex in range(validSmoothingIterations):
        validWeightSum = convolve2d(validWeightSum, validKernel, mode='same')
        validWeightSum[iterMask] = 0.

    invalidWeightSum = convolve2d(fillMask, invalidKernel, mode='same')

    ny, nx = fillMask.shape
    indices = numpy.indices((ny, nx))
    yIndices = indices[0].ravel()
    xIndices = indices[1].ravel()

    fillIndices = numpy.nonzero(fillMask.ravel())[0]
    invalidWeightSum = invalidWeightSum[fillMask]
    weightSum = validWeightSum[fillMask] + invalidWeightSum
    validWeightSum = validWeightSum[valid]

    fillInverseMap = -1*numpy.ones((ny, nx), int)
    fillInverseMap[fillMask] = numpy.arange(fillCount)

    matrix = lil_matrix((fillCount, fillCount))
    for index, fillIndex in enumerate(fillIndices):
        xc = xIndices[fillIndex]
        yc = yIndices[fillIndex]
        xMin = max(0, xc-invalidKernelSize)
        xMax = min(nx, xc+invalidKernelSize+1)
        yMin = max(0, yc-invalidKernelSize)
        yMax = min(ny, yc+invalidKernelSize+1)
        kxMin = xMin-xc+invalidKernelSize
        kxMax = xMax-xc+invalidKernelSize
        kyMin = yMin-yc+invalidKernelSize
        kyMax = yMax-yc+invalidKernelSize
        otherIndices = fillInverseMap[yMin:yMax, xMin:xMax]
        weights = invalidKernel[kyMin:kyMax, kxMin:kxMax]/weightSum[index]
        mask = otherIndices >= 0
        otherIndices = otherIndices[mask]
        weights = weights[mask]
        matrix[index, otherIndices] = -weights
        # add ones along the diagonal
        matrix[index, index] = 1 + matrix[index, index]

    matrix = matrix.tocsr()
    _save_matrix_and_kernel(outFileName, matrix, validKernel, valid,
                            fillMask, weightSum, validWeightSum)


def _extrap_basin(ds, fieldName, basinName, matrixFileTemplate, parallelTasks,
                  basinMask, validSmoothingIterations,
                  replaceValidWithSmoothed):

    nz = ds.sizes['z']

    field3D = ds[fieldName].values

    origShape = field3D.shape

    if 'time' not in ds.dims:
        field3D = field3D.reshape((1, origShape[0], origShape[1],
                                   origShape[2]))

    print('  Extrapolating {} in {}...'.format(fieldName, basinName))
    widgets = ['  ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nz).start()

    pool = Pool(parallelTasks)
    zIndices = range(nz)
    partial_func = partial(_extrap_basin_level, field3D, matrixFileTemplate,
                           basinMask, replaceValidWithSmoothed,
                           validSmoothingIterations)

    for zIndex, fieldSlice in enumerate(pool.imap(partial_func, zIndices)):
        field3D[:, zIndex, :, :] = fieldSlice
        bar.update(zIndex+1)

    bar.finish()

    if 'time' not in ds.dims:
        field3D = field3D.reshape(origShape)

    dims = ds[fieldName].dims
    attrs = ds[fieldName].attrs
    ds[fieldName] = (dims, field3D)
    ds[fieldName].attrs = attrs


def _extrap_basin_level(field3D, matrixFileTemplate, basinMask,
                        replaceValidWithSmoothed, validSmoothingIterations,
                        zIndex):

    matrix, validKernel, valid, fillMask, weightSum, validWeightSum = \
        _load_matrix_and_kernel(matrixFileTemplate.format(zIndex))

    nt, nz, ny, nx = field3D.shape

    nanMask = numpy.logical_not(numpy.logical_or(valid, fillMask))

    outField = field3D[:, zIndex, :, :]

    basinFillMask = numpy.logical_and(fillMask, basinMask)
    # no point in doing extrapolation if there are no fill points in the limits
    # of the basin
    basinFillCount = numpy.count_nonzero(basinFillMask)
    validCount = numpy.count_nonzero(valid)

    iterMask = numpy.logical_not(numpy.logical_or(valid, fillMask))

    for tIndex in range(nt):
        fieldSlice = outField[tIndex, :, :]
        if basinFillCount > 0 and validCount > 0:
            fieldExtrap = fieldSlice.copy()
            fieldExtrap[numpy.logical_not(valid)] = 0.
            fieldExtrap[numpy.isnan(fieldExtrap)] = 0.
            for iterIndex in range(validSmoothingIterations):
                fieldExtrap = convolve2d(fieldExtrap, validKernel, mode='same')
                fieldExtrap[iterMask] = 0.

            if replaceValidWithSmoothed:
                fieldSlice[valid] = fieldExtrap[valid]/validWeightSum
            rhs = fieldExtrap[fillMask]/weightSum

            fieldFill = spsolve(matrix, rhs)
            fieldSlice[fillMask] = fieldFill
        else:
            fieldSlice[fillMask] = numpy.nan

        fieldSlice[nanMask] = numpy.nan

        outField[tIndex, :, :] = fieldSlice

    return outField


def _add_basin_field(dsIn, dsOut, fieldName, basinMask):

    nz = dsIn.sizes['z']

    fieldIn = dsIn[fieldName].values
    fieldOut = dsOut[fieldName].values

    origShape = fieldIn.shape

    if 'time' in dsIn.dims:
        nt = dsIn.sizes['time']
    else:
        nt = 1
        fieldIn = fieldIn.reshape((1, origShape[0], origShape[1],
                                   origShape[2]))
        fieldOut = fieldOut.reshape((1, origShape[0], origShape[1],
                                     origShape[2]))

    for zIndex in range(nz):
        for tIndex in range(nt):
            fieldSliceIn = fieldIn[tIndex, zIndex, :, :]
            fieldSliceOut = fieldOut[tIndex, zIndex, :, :]
            fieldSliceOut[basinMask] = fieldSliceIn[basinMask]
            fieldOut[tIndex, zIndex, :, :] = fieldSliceOut

    if 'time' not in dsIn.dims:
        fieldOut = fieldOut.reshape(origShape)

    dsOut[fieldName] = (dsIn[fieldName].dims, fieldOut)
    dsOut[fieldName].attrs = dsIn[fieldName].attrs


def _save_matrix_and_kernel(fileName, matrix, kernel, valid, fillMask,
                            weightSum, validWeightSum):
    numpy.savez(fileName, data=matrix.data, indices=matrix.indices,
                indptr=matrix.indptr, shape=matrix.shape, kernel=kernel,
                valid=valid, fillMask=fillMask, weightSum=weightSum,
                validWeightSum=validWeightSum)


def _load_matrix_and_kernel(fileName):
    loader = numpy.load(fileName)
    kernel = loader['kernel']
    valid = loader['valid']
    fillMask = loader['fillMask']
    weightSum = loader['weightSum']
    validWeightSum = loader['validWeightSum']
    matrix = csr_matrix((loader['data'], loader['indices'], loader['indptr']),
                        shape=loader['shape'])
    return matrix, kernel, valid, fillMask, weightSum, validWeightSum
