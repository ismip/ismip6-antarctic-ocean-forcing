from ismip6_ocean_forcing.remap.grid import ProjectionGridDescriptor, \
    LatLonGridDescriptor
from ismip6_ocean_forcing.remap.polar import \
    get_antarctic_stereographic_projection

import xarray


def get_antarctic_descriptor(fileName):
    dsIn = xarray.open_dataset(fileName)
    x = dsIn.x.values
    y = dsIn.y.values
    dx = int((x[1]-x[0])/1000.)
    Lx = int((x[-1] - x[0])/1000.)
    Ly = int((y[-1] - y[0])/1000.)

    meshName = '{}x{}km_{}km_Antarctic_stereo'.format(Lx, Ly, dx)

    projection = get_antarctic_stereographic_projection()

    descriptor = ProjectionGridDescriptor.create(projection, x, y, meshName)

    return descriptor


def get_lat_lon_descriptor(fileName):

    descriptor = LatLonGridDescriptor.read(fileName, latVarName='lat',
                                           lonVarName='lon')

    return descriptor
