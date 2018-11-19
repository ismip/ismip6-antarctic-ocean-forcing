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

    lat = numpy.maximum(dsTemp.lat.values, -80.)
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
            t_freezing = gsw.freezing.t_freezing(SA, pressure[mask], 0.)

            thermalForcingLocal = numpy.nan*numpy.ones(temp.shape)
            thermalForcingLocal[mask] = temp[mask] - t_freezing
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


def potential_to_in_situ_temperature(dsPotTemp, dsSalin):
    z = dsPotTemp.z.values
    lat = numpy.maximum(dsPotTemp.lat.values, -80.)
    lon = dsPotTemp.lon.values

    nz = len(z)
    ny, nx = lat.shape

    if 'time' in dsPotTemp.dims:
        nt = dsPotTemp.sizes['time']
        T = numpy.nan*numpy.ones((nt, nz, ny, nx))
        for zIndex in range(nz):
            pressure = gsw.p_from_z(z[zIndex], lat)
            for tIndex in range(nt):
                pt = dsPotTemp.temperature[tIndex, zIndex, :, :].values
                salin = dsSalin.salinity[tIndex, zIndex, :, :].values
                mask = numpy.logical_and(numpy.isfinite(pt), numpy.isfinite(salin))
                SA = gsw.SA_from_SP(salin[mask], pressure[mask], lon[mask],
                                    lat[mask])
                TSlice = T[tIndex, zIndex, :, :]
                CT = gsw.CT_from_pt(SA, pt[mask])
                TSlice[mask] = gsw.t_from_CT(SA, CT, pressure[mask])
                T[tIndex, zIndex, :, :] = TSlice
    else:
        T = numpy.nan*numpy.ones((nz, ny, nx))
        for zIndex in range(nz):
            pressure = gsw.p_from_z(z[zIndex], lat)
            pt = dsPotTemp.temperature[zIndex, :, :].values
            salin = dsSalin.salinity[zIndex, :, :].values
            mask = numpy.logical_and(numpy.isfinite(pt), numpy.isfinite(salin))
            SA = gsw.SA_from_SP(salin[mask], pressure[mask], lon[mask],
                                lat[mask])
            TSlice = T[zIndex, :, :]
            CT = gsw.CT_from_pt(SA, pt[mask])
            TSlice[mask] = gsw.t_from_CT(SA, CT, pressure[mask])
            T[zIndex, :, :] = TSlice

    dsTemp = dsPotTemp.drop('temperature')
    dsTemp['temperature'] = (dsPotTemp.temperature.dims, T)
    dsTemp['temperature'].attrs = dsPotTemp.temperature.attrs

    return dsTemp
