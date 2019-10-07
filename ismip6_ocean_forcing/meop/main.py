import os
import xarray
import numpy
import glob
import gsw
import progressbar
import subprocess
from pyremap.polar import get_antarctic_stereographic_projection

from ismip6_ocean_forcing.io import download_files
from ismip6_ocean_forcing.remap.res import get_res, get_horiz_res


def process_meop(config):
    '''
    Download MEOP temperature and salinity, and bin them on the ISMIP6 grid
    '''

    try:
        os.makedirs('meop')
    except OSError:
        pass

    print('Processing Marine Mammals Exploring the Oceans Pole to Pole '
          '(MEOP)...')

    if not os.path.exists('meop/MEOP-CTD_2018-04-10'):
        baseURL = 'https://www.seanoe.org/data/00343/45461/data'
        fileNames = ['58202.zip']
        download_files(fileNames, baseURL, 'meop')

        print('  Decompressing MEOP data...')
        args = ['unzip', 'meop/58202.zip', '-d', 'meop/']
        returncode = subprocess.call(args)
        if returncode not in [0, 1, 2]:
            raise subprocess.CalledProcessError(returncode, args)

        print('     Done.')

    _bin_meop(config, 'TEMP', 'temperature')
    _bin_meop(config, 'PSAL', 'salinity')
    print('Done.')


def _bin_meop(config, inVarName, outVarName):
    res = get_res(config)
    outFileName = 'meop/meop_{}_{}.nc'.format(outVarName, res)
    if os.path.exists(outFileName):
        return

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

    fileList = sorted(glob.glob('meop/MEOP-CTD_2018-04-10/*/DATA_ncARGO/*.nc'))
    print('  Binning MEOP {} profiles...'.format(outVarName))

    widgets = ['  ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=len(fileList)).start()

    for index, fileName in enumerate(fileList):
        dsProfile = xarray.open_dataset(fileName)
        lat = dsProfile.LATITUDE.values
        lon = dsProfile.LONGITUDE.values
        inField = dsProfile['{}_ADJUSTED'.format(inVarName)].values
        quality = dsProfile['{}_ADJUSTED_QC'.format(inVarName)].values
        if attrs is None:
            attrs = dsProfile[inVarName].attrs
        x, y = proj(lon, lat)
        pressure = dsProfile.PRES.values

        lat = numpy.maximum(lat, -75.)
        for profile in range(pressure.shape[0]):
            xBin = int((x[profile]-xMin)/dx)
            yBin = int((y[profile]-yMin)/dx)
            if xBin < 0 or xBin >= nx:
                continue
            if yBin < 0 or yBin >= ny:
                continue
            for level in range(pressure.shape[1]):
                if quality[profile, level] != b'1':
                    continue
                press = pressure[profile, level]
                if numpy.isnan(press):
                    continue
                depth = gsw.z_from_p(pressure[profile, level], lat[profile])
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
    ds.to_netcdf(outFileName)
