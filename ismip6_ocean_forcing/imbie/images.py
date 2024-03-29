import shapefile
import xarray
import shapely.geometry
import shapely.ops
from descartes import PolygonPatch
import os.path
import matplotlib.pyplot as plt


def write_basin_images(res, inFileName, basins):

    if os.path.exists('imbie/basins_{}/'.format(res)):
        return

    basinFileName = 'imbie/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'

    print('  Loading basin geometry...')
    reader = shapefile.Reader(basinFileName)
    fields = reader.fields[1:]
    field_names = [field[0] for field in fields]
    inBasinData = {}
    for sr in reader.shapeRecords():
        atr = dict(zip(field_names, sr.record))
        if atr['Subregion'] == '':
            continue
        atr['name'] = atr['Subregion']
        geom = sr.shape.__geo_interface__
        inBasinData[atr['name']] = dict(type="Feature",
                                        geometry=geom, properties=atr)

    outBasinData = []
    for outBasinName in basins:
        inBasins = basins[outBasinName]
        if len(inBasins) == 1:
            inBasinName = inBasins[0]
            outBasinData.append(inBasinData[inBasinName])
        else:
            featuresToCombine = []
            for inBasinName in inBasins:
                featuresToCombine.append(inBasinData[inBasinName])
            feature = _combine_features(featuresToCombine, outBasinName)
            outBasinData.append(feature)

    ds = xarray.open_dataset(inFileName)
    nx = ds.sizes['x']
    ny = ds.sizes['y']

    dx = (ds.x[1] - ds.x[0]).values

    try:
        os.makedirs('imbie/basins_{}'.format(res))
    except OSError:
        pass

    print('  Writing basin images...')
    for feature in outBasinData:
        name = feature['properties']['name']
        basinShape = shapely.geometry.shape(feature['geometry'])
        print('    {}'.format(name))
        _write_basin_image(res, basinShape, name, nx, ny, dx)


def _combine_features(featuresToCombine, newName):
    featureShapes = []
    featureNames = []
    for feature in featuresToCombine:
        featureShapes.append(shapely.geometry.shape(feature['geometry']))
        featureNames.append(feature['properties']['name'])

    combinedShape = shapely.ops.unary_union(featureShapes)

    feature = {}
    feature['properties'] = {}
    feature['properties']['name'] = newName
    feature['geometry'] = shapely.geometry.mapping(combinedShape)

    if feature['geometry']['type'] == 'GeometryCollection':
        raise ValueError(
            "Error: combined geometry is of type GeometryCollection.\n"
            "       Most likely cause is that multiple feature types "
            "(regions, \n"
            "       points and transects) are being cobined.")

    return feature


def _write_basin_image(res, basinShape, name, nx, ny, dx):
    my_dpi = 600
    fig = plt.figure(figsize=(nx/float(my_dpi), ny/float(my_dpi)), dpi=my_dpi)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    plt.axis('off')

    color = 'black'

    if basinShape.geom_type == 'Polygon':
        ax.add_patch(PolygonPatch(
            basinShape.__geo_interface__, fc=color, ec=color, linewidth=0))
    elif basinShape.geom_type == 'MultiPolygon':
        for poly in basinShape:
            ax.add_patch(PolygonPatch(
                poly.__geo_interface__, fc=color, ec=color, linewidth=0))

    plt.xlim([-dx*0.5*(nx+1), dx*0.5*(nx+1)])
    plt.ylim([-dx*0.5*(ny+1), dx*0.5*(ny+1)])
    fig.canvas.draw()
    plt.savefig('imbie/basins_{}/{}.png'.format(res, name), dpi=my_dpi)
    plt.close()
