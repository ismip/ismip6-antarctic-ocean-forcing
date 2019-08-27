Download and Process CESM2 thetao and so
========================================

These data are available on the ISMIP6 lftp server.  The files are in the
subdirectory `/GSFC/Kate`.  The files needed are:

* `b.e21.BHIST.f09_g17.CMIP6-historical.011.pop.h.SALT.195001-199912.nc`
* `b.e21.BHIST.f09_g17.CMIP6-historical.011.pop.h.SALT.200001-201412.nc`
* `b.e21.BHIST.f09_g17.CMIP6-historical.011.pop.h.TEMP.195001-199912.nc`
* `b.e21.BHIST.f09_g17.CMIP6-historical.011.pop.h.TEMP.200001-201412.nc`
* `b.e21.BSSP585cmip6.f09_g17.CMIP6-SSP5-8.5.002.pop.h.SALT.201501-206412.nc`
* `b.e21.BSSP585cmip6.f09_g17.CMIP6-SSP5-8.5.002.pop.h.SALT.206501-210012.nc`
* `b.e21.BSSP585cmip6.f09_g17.CMIP6-SSP5-8.5.002.pop.h.TEMP.201501-206412.nc`
* `b.e21.BSSP585cmip6.f09_g17.CMIP6-SSP5-8.5.002.pop.h.TEMP.206501-210012.nc`

Finally, run `process_cesm2_e21_f09_g17.py -o /path/to/cesm2/files/` with the
path where the CESM2 data sets are stored.

