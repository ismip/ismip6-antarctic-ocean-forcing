import xarray
import numpy
import os
import progressbar
import gsw
import gsw.freezing


def compute_thermal_forcing(temperatureFileName, salinityFileName, outFileName,
                            timeIndices=None):

    if os.path.exists(outFileName):
        return

    dsTemp = xarray.open_dataset(temperatureFileName)
    dsSalin = xarray.open_dataset(salinityFileName)
    if 'time' in dsTemp.dims and timeIndices is not None:
        dsTemp = dsTemp.isel(time=timeIndices)
        dsSalin = dsSalin.isel(time=timeIndices)

    nx = dsTemp.sizes['x']
    ny = dsTemp.sizes['y']
    nz = dsTemp.sizes['z']

    if 'time' in dsTemp.dims:
        nt = dsTemp.sizes['time']
    else:
        nt = 1

    thermalForcing = numpy.zeros((nt, nz, ny, nx))

    print('  Computing thermal forcing...')
    widgets = ['  z=1/{}: '.format(nz),
               progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nz*nt).start()

    lat = dsTemp.lat.values
    lon = dsTemp.lon.values
    for zIndex in range(nz):
        pressure = gsw.p_from_z(dsTemp.z[zIndex].values, lat)
        for tIndex in range(nt):

            if 'time' in dsTemp.dims:
                temp = dsTemp.temperature[tIndex, zIndex, :, :].values
                salin = dsSalin.salinity[tIndex, zIndex, :, :].values
            else:
                temp = dsTemp.temperature[zIndex, :, :].values
                salin = dsSalin.salinity[zIndex, :, :].values
            mask = numpy.isfinite(temp)
            SA = gsw.SA_from_SP(salin[mask], pressure[mask], lon[mask],
                                lat[mask])
            CT = gsw.CT_from_pt(SA, temp[mask])
            CT_freezing = gsw.freezing.CT_freezing(SA, pressure[mask], 0.)

            thermalForcingLocal = numpy.nan*numpy.ones(temp.shape)
            thermalForcingLocal[mask] = CT - CT_freezing
            thermalForcing[tIndex, zIndex, :, :] = thermalForcingLocal

            widgets[0] = '  z={}/{}: '.format(zIndex+1, nz)
            bar.update(tIndex + nt*zIndex)

    bar.finish()

    if 'time' not in dsTemp.dims:
        thermalForcing = thermalForcing.reshape((nz, ny, nx))

    fieldName = 'thermal_forcing'
    ds = xarray.Dataset()
    for var in ['x', 'y', 'z', 'z_bnds']:
        ds[var] = dsTemp[var]
    if 'time' in dsTemp.dims:
        ds['time'] = dsTemp['time']

    ds[fieldName] = (dsTemp.temperature.dims, thermalForcing)
    ds[fieldName].attrs['units'] = "degrees_celsius"
    ds[fieldName].attrs['long_name'] = "thermal forcing, the difference " \
        "between the temperature and the local freezing point"

    ds.to_netcdf(outFileName)
