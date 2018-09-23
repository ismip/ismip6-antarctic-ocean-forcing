import pyproj


def get_antarctic_stereographic_projection():
    """
    Get a projection for an Antarctic steregraphic grid
    """

    projection = pyproj.Proj('+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 '
                             '+k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84')

    return projection


def to_polar(points):

    projection = get_antarctic_stereographic_projection()
    latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')

    x, y = pyproj.transform(latLonProjection, projection, points[:, 0],
                            points[:, 1], radians=False)
    points[:, 0] = x
    points[:, 1] = y
    return points


def from_polar(points):

    projection = get_antarctic_stereographic_projection()
    latLonProjection = pyproj.Proj(proj='latlong', datum='WGS84')

    lon, lat = pyproj.transform(projection, latLonProjection, points[:, 0],
                                points[:, 1], radians=False)
    points[:, 0] = lon
    points[:, 1] = lat
    return points
