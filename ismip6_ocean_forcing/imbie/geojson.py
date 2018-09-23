import json
import shapely.geometry
import shapely.ops
from collections import OrderedDict
import os

from ismip6_ocean_forcing.io import download_files


def download_imbie():
    '''
    Download the geojson files that define the IMBIE basins
    '''

    geojsonURL = 'https://raw.githubusercontent.com/MPAS-Dev/' \
                 'geometric_features/master/'

    fileNames = [
            'landice/region/Antarctica_IMBIE{}/region.geojson'.format(basin)
            for basin in range(1, 28)]

    download_files(fileNames, geojsonURL, 'imbie')


def combine_imbie(basins):

    outFileName = 'imbie/AntarcticBasins.geojson'
    if os.path.exists(outFileName):
        return

    print('  Combining IMBIE basins...')
    template = 'imbie/landice/region/Antarctica_IMBIE{}/region.geojson'

    features = {'features': []}
    for basinList in basins:
        if isinstance(basinList, list):
            fileNames = [template.format(basin) for basin in basinList]

            basinName = 'Antarctica_IMBIE{}'.format(
                    '_'.join(['{}'.format(basin) for basin in basinList]))

            print('    {}'.format(basinName))
            featuresToCombine = {'features': []}

            for fileName in fileNames:
                _merge_features(fileName, featuresToCombine)

            _combine_features(featuresToCombine, basinName, features)

        else:
            fileName = template.format(basinList)
            print('    Antarctica_IMBIE{}'.format(basinList))
            _merge_features(fileName, features)

    _write_geojson(features, outFileName)


def _merge_features(fileName, features):
    with open(fileName) as f:
        newFeatures = json.load(f)
        for feature in newFeatures['features']:
            features['features'].append(feature)


def _combine_features(featuresToCombine, newName, features):
    featureShapes = []
    authors = []
    featureNames = []
    for feature in featuresToCombine['features']:
        featureShapes.append(shapely.geometry.shape(feature['geometry']))
        authors.append(feature['properties']['author'])
        featureNames.append(feature['properties']['name'])

    combinedShape = shapely.ops.cascaded_union(featureShapes)

    feature = {}
    feature['properties'] = {}
    feature['properties']['name'] = newName
    feature['properties']['component'] = \
        featuresToCombine['features'][0]['properties']['component']
    feature['properties']['tags'] = ''
    feature['properties']['author'] = '; '.join(list(set(authors)))
    feature['properties']['constituents'] = '; '.join(list(set(featureNames)))
    feature['geometry'] = shapely.geometry.mapping(combinedShape)

    if feature['geometry']['type'] == 'GeometryCollection':
        raise ValueError(
            "Error: combined geometry is of type GeometryCollection.\n"
            "       Most likely cause is that multiple feature types "
            "(regions, \n"
            "       points and transects) are being cobined.")

    features['features'].append(feature)


def _write_geojson(features, fileName):

    features['type'] = 'FeatureCollection'

    for index in range(len(features['features'])):
        features['features'][index] = \
            _check_feature(features['features'][index])

    # Make the feature an ordered dictionary so type comes before
    # features (easier to read)
    outFeatures = OrderedDict((('type', features['type']),))

    # Add the rest (except features)
    for key in sorted(features):
        if key not in ['type', 'features']:
            outFeatures[key] = features[key]

    # features go last for readability
    outFeatures['features'] = features['features']

    out_file = open(fileName, 'w')

    json.dump(outFeatures, out_file, indent=4)


def _check_feature(feature):  # {{{

    if 'name' not in feature['properties'].keys():
        raise ValueError("feature has no 'name' property.")

    if 'type' not in feature['geometry'].keys():
        raise ValueError("Feature {} has an issue with the geometry type."
                         "".format(feature['properties']['name']))

    featureType = feature['geometry']['type']

    # Determine object property value based on feature type.
    if featureType == "Polygon" or featureType == "MultiPolygon":
        objectType = "region"
    elif featureType == "LineString" or featureType == "MultiLineString":
        objectType = "transect"
    elif featureType == "Point" or featureType == "MultiPoint":
        objectType = "point"
    else:
        raise ValueError("Unsupported feature type {}".format(featureType))

    feature['properties']['object'] = objectType

    # Make the properties an ordered dictionary so they can be sorted
    outProperties = OrderedDict(
            (('name', feature['properties']['name']),
             ('object', feature['properties']['object'])))
    for key in sorted(feature['properties']):
        if key not in outProperties.keys():
            outProperties[key] = feature['properties'][key]

    # Make the geometry an ordered dictionary so they can keep it in the
    # desired order
    outGeometry = OrderedDict(
        (('type', feature['geometry']['type']),
         ('coordinates', feature['geometry']['coordinates'])))
    for key in sorted(feature['geometry']):
        if key not in outGeometry.keys():
            outGeometry[key] = feature['geometry'][key]

    # Make the feature an ordered dictionary so properties come before geometry
    # (easier to read)
    outFeature = OrderedDict((('type', 'Feature'),
                             ('properties', outProperties)))
    # Add the rest
    for key in sorted(feature):
        if key not in ['geometry', 'type', 'properties']:
            outFeature[key] = feature[key]

    # geometry goes last for readability
    outFeature['geometry'] = outGeometry

    return outFeature
