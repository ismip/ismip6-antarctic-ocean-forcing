import os
import xarray
import numpy
import glob
import progressbar
import zipfile
from pyremap.polar import get_antarctic_stereographic_projection

from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res
from ismip6_ocean_forcing.thermal_forcing.main import \
    potential_to_in_situ_temperature


def process_en4(config, startYear, endYear):
    '''
    Download EN4 temperature and salinity, and bin them on the ISMIP6 grid
    '''

    try:
        os.makedirs('en4/zips')
    except OSError:
        pass

    print('Processing UK Met Office EN4...')

    if not os.path.exists('en4/profiles'):
        baseURL = 'https://www.metoffice.gov.uk/hadobs/en4/data/en4-2-1/'
        fileNames = ['EN.4.2.1.profiles.g10.{}.zip'.format(year) for year in
                     range(startYear, endYear+1)]
        download_files(fileNames, baseURL, 'en4/zips')

        try:
            os.makedirs('en4/profiles')
        except OSError:
            pass

        print('  Decompressing EN4 data...')
        widgets = ['  ', progressbar.Percentage(), ' ',
                   progressbar.Bar(), ' ', progressbar.ETA()]
        bar = progressbar.ProgressBar(widgets=widgets,
                                      maxval=len(fileNames)).start()
        for index, fileName in enumerate(fileNames):
            with zipfile.ZipFile('en4/zips/{}'.format(fileName), 'r') as f:
                f.extractall('en4/profiles')
            bar.update(index+1)
        bar.finish()

    res = get_res(config)
    tempFileName = 'en4/en4_temperature_{}-{}_{}.nc'.format(startYear,
                                                            endYear, res)
    salinFileName = 'en4/en4_salinity_{}-{}_{}.nc'.format(startYear,
                                                          endYear, res)
    if not os.path.exists(salinFileName):
        dsSalinity = _bin_en4(config, 'PSAL', 'salinity', startYear, endYear)
        dsSalinity.to_netcdf(salinFileName)
        dsSalinity.close()

    if not os.path.exists(tempFileName):
        dsSalinity = xarray.open_dataset(salinFileName)
        dsPotTemp = _bin_en4(config, 'POTM', 'temperature', startYear, endYear)
        dsTemp = potential_to_in_situ_temperature(dsPotTemp, dsSalinity)
        dsTemp.to_netcdf(tempFileName)

    print('Done.')


def _bin_en4(config, inVarName, outVarName, startYear, endYear):
    hres = get_horiz_res(config)
    dz = config.getfloat('grid', 'dzExtrap')
    nz = config.getint('grid', 'nzExtrap')
    zOut = dz*numpy.arange(nz+1)
    z = 0.5*(zOut[0:-1] + zOut[1:])
    z_bnds = numpy.zeros((len(z), 2))
    z_bnds[:, 0] = zOut[0:-1]
    z_bnds[:, 1] = zOut[1:]

    ds = xarray.open_dataset('ismip6/{}_grid.nc'.format(hres))
    ds['z'] = (('z',), z)
    ds.z.attrs['units'] = 'meters'
    ds.z.attrs['bounds'] = 'z_bnds'
    ds.z.attrs['standard_name'] = 'depth'
    ds.z.attrs['positive'] = 'up'
    ds.z.attrs['axis'] = 'Z'

    ds['z_bnds'] = (('z', 'nbounds'), z_bnds)
    ds.z_bnds.attrs['comment'] = 'depth bounds'

    xMin = ds.x[0].values
    yMin = ds.y[0].values
    zMin = z[0]
    dx = ds.x[1].values - ds.x[0].values

    nx = ds.sizes['x']
    ny = ds.sizes['y']
    nz = ds.sizes['z']

    outField = numpy.zeros((nz, ny, nx))
    entryCount = numpy.zeros((nz, ny, nx), dtype=int)

    attrs = None

    proj = get_antarctic_stereographic_projection()

    X, Y = numpy.meshgrid(ds.x.values, ds.y.values)
    Lon, Lat = proj(X, Y, inverse=True)
    ds['lon'] = (('y', 'x'), Lon)
    ds.lon.attrs['units'] = 'degrees'
    ds['lat'] = (('y', 'x'), Lat)
    ds.lat.attrs['units'] = 'degrees'

    fileList = sorted(glob.glob('en4/profiles/*.nc'))
    print('  Binning EN4 {} profiles...'.format(outVarName))

    widgets = ['  ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=len(fileList)).start()

    for index, fileName in enumerate(fileList):
        dsProfile = xarray.open_dataset(fileName)
        lat = dsProfile.LATITUDE.values
        lon = dsProfile.LONGITUDE.values
        inField = dsProfile['{}_CORRECTED'.format(inVarName)].values
        quality = dsProfile['{}_CORRECTED_QC'.format(inVarName)].values
        if attrs is None:
            attrs = dsProfile['{}_CORRECTED'.format(inVarName)].attrs
        x, y = proj(lon, lat)
        depths = dsProfile.DEPH_CORRECTED.values

        lat = numpy.maximum(lat, -75.)
        for profile in range(depths.shape[0]):
            xBin = int((x[profile]-xMin)/dx)
            yBin = int((y[profile]-yMin)/dx)
            if xBin < 0 or xBin >= nx:
                continue
            if yBin < 0 or yBin >= ny:
                continue
            for level in range(depths.shape[1]):
                if quality[profile, level] != b'1':
                    continue
                depth = depths[profile, level]
                if numpy.isnan(depth):
                    continue
                zBin = int((depth-zMin)/dz)
                if zBin < 0 or zBin >= nz:
                    continue
                outField[zBin, yBin, xBin] += inField[profile, level]
                entryCount[zBin, yBin, xBin] += 1
        bar.update(index+1)
    bar.finish()
    mask = entryCount > 0
    outField[mask] /= entryCount[mask]
    outField[numpy.logical_not(mask)] = numpy.nan
    ds[outVarName] = (('z', 'y', 'x'), outField)
    for attr in ['units', 'long_name', 'comment']:
        ds[outVarName].attrs[attr] = attrs[attr]
    return ds
