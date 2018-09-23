import numpy
import json
import xarray
from shapely.geometry import shape
from descartes import PolygonPatch
import os.path
import matplotlib.pyplot as plt
from ismip6_ocean_forcing.remap.polar import to_polar


def write_basin_images(inFileName):

    if os.path.exists('imbie/basins/'):
        return

    basinFileName = 'imbie/AntarcticBasins.geojson'

    print('  Loading basin geometry...')
    with open(basinFileName) as f:
        basinData = json.load(f)

    print('  Converting basins to polar coordintes...')
    basinShapes = _make_polar_basins(basinData)

    ds = xarray.open_dataset(inFileName)
    nx = ds.sizes['x']
    ny = ds.sizes['y']

    dx = (ds.x[1] - ds.x[0]).values

    try:
        os.makedirs('imbie/basins')
    except OSError:
        pass

    print('  Writing basin images...')
    for index in range(len(basinShapes)):
        name = basinData['features'][index]['properties']['name']
        print('    {}'.format(name))
        _write_basin_image(basinShapes[index], name, nx, ny, dx)


def _make_polar_basins(basinData):
    basinShapes = []
    for feature in basinData['features']:
        print('    {}'.format(feature['properties']['name']))
        basinGeom = feature['geometry']
        coords = basinGeom['coordinates']
        newCoords = []
        if(basinGeom['type'] == 'Polygon'):
            for subpoly in coords:
                newCoords.append(to_polar(numpy.array(subpoly)).tolist())
        elif(basinGeom['type'] == 'MultiPolygon'):
            for poly in coords:
                newPoly = []
                for subpoly in poly:
                    newPoly.append(to_polar(numpy.array(subpoly)).tolist())
                newCoords.append(newPoly)
        basinGeom['coordinates'] = newCoords
        basinShape = shape(basinGeom)
        basinShapes.append(basinShape)
    return basinShapes


def _write_basin_image(basinShape, name, nx, ny, dx):
    my_dpi = 600
    fig = plt.figure(figsize=(nx/float(my_dpi), ny/float(my_dpi)), dpi=my_dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.axis('off')

    color = 'black'

    if basinShape.geom_type == 'Polygon':
        ax.add_patch(PolygonPatch(basinShape, fc=color, ec=color, linewidth=0))
    elif basinShape.geom_type == 'MultiPolygon':
        for poly in basinShape:
            ax.add_patch(PolygonPatch(poly, fc=color, ec=color, linewidth=0))

    plt.xlim([-dx*0.5*(nx+1), dx*0.5*(nx+1)])
    plt.ylim([-dx*0.5*(ny+1), dx*0.5*(ny+1)])
    fig.canvas.draw()
    plt.savefig('imbie/basins/{}.png'.format(name), dpi=my_dpi)
    plt.close()
