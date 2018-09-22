import os
import xarray
from ismip6_ocean_forcing.thermal_forcing.main import compute_thermal_forcing
from ismip6_ocean_forcing.remap.res import get_res


def compute_anomaly_and_to_woa(config):
    '''
    Combine data from various analysis runs into a single data set, comput its
    anomaly with respect to a present-day climatology, add the result to
    the WOA climatology, and compute the thermal forcing
    '''

    if not config.getboolean('combine', 'combine'):
        return

    modelName = config.get('model', 'name')
    print('Computing anomalies of {} data and adding them to WOA...'.format(
            modelName))

    _combine_model_output(config)
    _compute_climatology(config)
    _compute_anomaly(config)
    _add_anomaly_to_woa(config)
    _compute_thermal_driving(config)
    print('  Done.')


def _combine_model_output(config):

    res = get_res(config)
    modelName = config.get('model', 'name')
    folders = config.get('combine', 'folders')
    combineDim = config.get('combine', 'dim')
    baseFolder = modelName.lower()
    outFolder = '{}/{}'.format(baseFolder, config.get('combine', 'outFolder'))
    folders = ['{}/{}'.format(baseFolder, folder.strip()) for folder in
               folders.split(',')]

    try:
        os.makedirs('outFolder')
    except OSError:
        pass

    print('  Combining model results into a single data set...')

    for fieldName in ['temperature', 'salinity']:
        fileNames = ['{}/{}_{}_{}_extrap_vert.nc'.format(
                folder, modelName, fieldName, res) for folder in folders]
        outFileName = '{}/{}_{}_{}.nc'.format(
                outFolder, modelName, fieldName, res)

        if os.path.exists(outFileName):
            print('    {} already exists.'.format(outFileName))
            continue

        print('    {}'.format(outFileName))
        ds = xarray.open_mfdataset(fileNames, concat_dim=combineDim)
        ds.to_netcdf(outFileName)


def _compute_climatology(config):
    res = get_res(config)
    modelName = config.get('model', 'name')
    firstTIndex = config.getint('climatology', 'firstTIndex')
    lastTIndex = config.getint('climatology', 'lastTIndex')
    baseFolder = modelName.lower()
    inFolder = '{}/{}'.format(baseFolder, config.get('combine', 'outFolder'))
    outFolder = '{}/{}'.format(baseFolder, config.get('climatology', 'folder'))

    print('  Computing present-day climatology...')
    for fieldName in ['temperature', 'salinity']:
        inFileName = '{}/{}_{}_{}.nc'.format(
                inFolder, modelName, fieldName, res)
        outFileName = '{}/{}_{}_{}.nc'.format(
                outFolder, modelName, fieldName, res)

        if os.path.exists(outFileName):
            print('    {} already exists.'.format(outFileName))
            continue

        print('    {}'.format(outFileName))
        ds = xarray.open_dataset(inFileName)
        ds = ds.isel(time=slice(firstTIndex, lastTIndex+1))
        ds = ds.mean(dim='time')
        ds.to_netcdf(outFileName)


def _compute_anomaly(config):
    res = get_res(config)
    modelName = config.get('model', 'name')
    baseFolder = modelName.lower()
    inFolder = '{}/{}'.format(baseFolder, config.get('combine', 'outFolder'))
    climFolder = '{}/{}'.format(baseFolder,
                                config.get('climatology', 'folder'))
    outFolder = '{}/{}'.format(baseFolder, config.get('anomaly', 'folder'))

    print('  Computing anomaly from present-day...')
    for fieldName in ['temperature', 'salinity']:
        inFileName = '{}/{}_{}_{}.nc'.format(
                inFolder, modelName, fieldName, res)
        climFileName = '{}/{}_{}_{}.nc'.format(
                climFolder, modelName, fieldName, res)
        outFileName = '{}/{}_{}_{}.nc'.format(
                outFolder, modelName, fieldName, res)

        if os.path.exists(outFileName):
            print('    {} already exists.'.format(outFileName))
            continue

        print('    {}'.format(outFileName))
        ds = xarray.open_dataset(inFileName)
        dsClim = xarray.open_dataset(climFileName)
        attrs = ds[fieldName].attrs
        ds[fieldName] = ds[fieldName] - dsClim[fieldName]
        ds[fieldName].attrs = attrs
        ds.to_netcdf(outFileName)


def _add_anomaly_to_woa(config):
    res = get_res(config)
    modelName = config.get('model', 'name')
    baseFolder = modelName.lower()
    inFolder = '{}/{}'.format(baseFolder, config.get('anomaly', 'folder'))
    outFolder = '{}/{}'.format(baseFolder, config.get('anomaly', 'woaFolder'))

    print('  Adding WOA climatology to the anomaly...')
    for fieldName in ['temperature', 'salinity']:
        woaFileName = \
            'woa/woa_{}_1995-2012_{}_extrap_vert.nc'.format(fieldName, res)
        inFileName = '{}/{}_{}_{}.nc'.format(
                inFolder, modelName, fieldName, res)
        outFileName = '{}/{}_{}_{}.nc'.format(
                outFolder, modelName, fieldName, res)

        if os.path.exists(outFileName):
            print('    {} already exists.'.format(outFileName))
            continue

        print('    {}'.format(outFileName))
        ds = xarray.open_dataset(inFileName)
        dsWOA = xarray.open_dataset(woaFileName)
        attrs = ds[fieldName].attrs
        ds[fieldName] = dsWOA[fieldName] + ds[fieldName]
        ds[fieldName].attrs = attrs
        ds.to_netcdf(outFileName)


def _compute_thermal_driving(config):
    res = get_res(config)
    modelName = config.get('model', 'name')
    subfolder = config.get('anomaly', 'woaFolder')
    modelFolder = '{}/{}'.format(modelName.lower(), subfolder)

    tempFileName = '{}/{}_temperature_{}.nc'.format(
            modelFolder, modelName, res)
    salinFileName = '{}/{}_salinity_{}.nc'.format(
            modelFolder, modelName, res)
    outFileName = '{}/{}_thermal_forcing_{}.nc'.format(
            modelFolder, modelName, res)
    compute_thermal_forcing(tempFileName, salinFileName, outFileName)
