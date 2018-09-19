import numpy
import xarray


def weights_and_indices(xInCenter=None, xOutCenter=None,
                        xInBounds=None, xOutBounds=None,
                        xDim='xOut', slabDim='slab'):
    '''
    Get the weights and indices for performing 1D conservative interpolation.
    Either center or bounds locations should be specified for each grid.

    The resuling weights and indices are given in slabs.  Interpolation is
    performed by summing over slabs the product of the given weights with
    the field indexed at the given indices::

        fieldOut = (field[inIndices]*weights).sum(dim=slabDim)

    Parameters
    ----------
    xInCenter : 1D float array, optional
        The locations of centers on the input grid, from which bounds will be
        interpolated or extrapolated

    xOutCenter : 1D float array, optional
        The locations of centers on the output grid, from which bounds will be
        interpolated or extrapolated

    xInBounds : 1D float array, optional
        The locations of bounds for each cell on the input grid

    xOutBounds : 1D float array, optional
        The locations of bounds for each cell on the output grid

    xDim : str, optional
        The name of the new dimension after interpolation

    slabDim : str, optional
        The name of the "slab" dimension to sum over when performing
        the interpolation

    Returns
    -------
    weightSlabs : xarray.DataArray
        Array of weights in "slabs", each the same size as xOutCenter.
    inIndexSlabs : xarray.DataArray
        Array of indices into the input array in "slabs", each the same size
        as xOutCenter.
    '''
    if xInBounds is None:
        xInBounds = numpy.zeros(len(xInCenter)+1)
        xInBounds[0] = 1.5*xInCenter[0] - 0.5*xInCenter[1]
        xInBounds[1:-1] = 0.5*(xInCenter[0:-1] + xInCenter[1:])
        xInBounds[-1] = 1.5*xInCenter[-1] - 0.5*xInCenter[-2]

    if xOutBounds is None:
        xOutBounds = numpy.zeros(len(xOutCenter)+1)
        xOutBounds[0] = 1.5*xOutCenter[0] - 0.5*xOutCenter[1]
        xOutBounds[1:-1] = 0.5*(xOutCenter[0:-1] + xOutCenter[1:])
        xOutBounds[-1] = 1.5*xOutCenter[-1] - 0.5*xOutCenter[-2]

    intersections = []
    nIn = len(xInBounds)-1
    nOut = len(xOutBounds)-1
    count = numpy.zeros(nOut, int)
    for iOut in range(nOut):
        dxOut = xOutBounds[iOut+1] - xOutBounds[iOut]
        for iIn in range(nIn):
            if dxOut > 0:
                xMin = max(xInBounds[iIn], xOutBounds[iOut])
                xMax = min(xInBounds[iIn+1], xOutBounds[iOut+1])
            else:
                xMin = min(xInBounds[iIn], xOutBounds[iOut])
                xMax = max(xInBounds[iIn+1], xOutBounds[iOut+1])
            weight = (xMax - xMin)/dxOut

            if weight <= 0.:
                continue
            intersections.append((iIn, iOut, weight))
            count[iOut] += 1

    maxCount = numpy.amax(count)
    weightSlabs = numpy.zeros((maxCount, nOut))
    inIndexSlabs = numpy.zeros((maxCount, nOut), int)
    slabIndex = numpy.zeros(nOut, int)
    for iIn, iOut, weight in intersections:
        sIndex = slabIndex[iOut]
        weightSlabs[sIndex, iOut] = weight
        inIndexSlabs[sIndex, iOut] = iIn
        slabIndex[iOut] += 1

    weightSlabs = xarray.DataArray(weightSlabs, dims=(slabDim, xDim))
    inIndexSlabs = xarray.DataArray(inIndexSlabs, dims=(slabDim, xDim))

    return weightSlabs, inIndexSlabs


def interp1d(field, weights, inIndices, normalizationThreshold=None,
             slabDim='slab'):
    '''
    1D interpolation of a field with weights and indices produced by
    ``weights_and_indices()``

    Parameters
    ----------
    field : xarray.DataArray
        A 1D field to interpolate

    weights : xarray.DataArray
        Array of weights in "slabs"

    inIndices : xarray.DataArray
        Array of indices into ``field`` in "slabs"

    normalizationThreshold : float, optional
        A threshold indicating that the results should be normalized and
        locations with weights summing to less than the threshold should be
        masked.  By default, regions where the weights sum to zero are masked
        and no normalization is performed

    slabDim : str, optional
        The name of the slab dim over which to perform the interpolation sum

    Returns
    -------
    result : xarray.DataArray
        ``field`` after interpolation
    '''

    result = (field[inIndices]*weights).sum(dim=slabDim)

    mask = field.notnull().astype(float)
    outMask = (mask[inIndices]*weights).sum(dim=slabDim)

    result = _normalize(result, outMask, normalizationThreshold)

    return result


def interp2d(field, xWeights, xIndices, yWeights, yIndices,
             normalizationThreshold=None, slabDim='slab'):
    '''
    2D interpolation performed as a sequence of 2 1D interpolations (on a
    "plaid" grid) of a field with weights and indices for both x and y
    produced by ``weights_and_indices()``

    Parameters
    ----------
    field : xarray.DataArray
        A 2D field to interpolate

    xWeights, yWeights : xarray.DataArray
        Array of weights in "slabs"

    xIndices, yIndices : xarray.DataArray
        Array of x and y indices into ``field`` in "slabs"

    normalizationThreshold : float, optional
        A threshold indicating that the results should be normalized and
        locations with weights summing to less than the threshold should be
        masked.  By default, regions where the weights sum to zero are masked
        and no normalization is performed

    slabDim : str, optional
        The name of the slab dim over which to perform the interpolation sum

    Returns
    -------
    result : xarray.DataArray
        ``field`` after interpolation
    '''
    temp = (field[:, xIndices]*xWeights).sum(dim=slabDim)
    result = (temp[yIndices, :]*yWeights).sum(dim=slabDim)

    mask = field.notnull().astype(float)
    temp = (mask[:, xIndices]*xWeights).sum(dim=slabDim)
    outMask = (temp[yIndices, :]*yWeights).sum(dim=slabDim)

    result = _normalize(result, outMask, normalizationThreshold)

    return result


def interp_depth(field, weights, inIndices, normalizationThreshold=None):
    '''
    Like ``interp1d`` but for the first dimension of a 3D array

    Parameters
    ----------
    field : xarray.DataArray
        A 1D field to interpolate

    weights : xarray.DataArray
        Array of weights in "slabs"

    inIndices : xarray.DataArray
        Array of indices into ``field`` in "slabs"

    normalizationThreshold : float, optional
        A threshold indicating that the results should be normalized and
        locations with weights summing to less than the threshold should be
        masked.  By default, regions where the weights sum to zero are masked
        and no normalization is performed

    Returns
    -------
    result : xarray.DataArray
        ``field`` after interpolation
    '''

    nzIn, ny, nx = field.values.shape
    nSlabs, nzOut = weights.values.shape

    outDims = (weights.dims[1], field.dims[1], field.dims[2])
    result = xarray.DataArray(numpy.zeros((nzOut, ny, nx)), dims=outDims)
    mask = xarray.DataArray(numpy.zeros((nzOut, ny, nx)), dims=outDims)
    for slab in range(nSlabs):
        fieldSlab = field[inIndices[slab, :], :, :]
        fieldSlab = fieldSlab.drop(field.dims[0])
        maskSlab = fieldSlab.notnull().astype(float)
        weightSlab = weights[slab, :]
        result += fieldSlab*weightSlab
        mask += maskSlab*weightSlab

    result = _normalize(result, mask, normalizationThreshold)

    return result


def _normalize(result, mask, normalizationThreshold):
    if normalizationThreshold is None:
        result = result.where(mask > 0.)
    else:
        # normalize by the mask and re-mask locations with too few valid
        # source points
        result = result.where(mask > normalizationThreshold)
        mask = mask.where(mask > normalizationThreshold, other=1.)
        result = result / mask

    return result
