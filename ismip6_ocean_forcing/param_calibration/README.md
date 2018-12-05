## Fortran scripts used to tune the melting parameterizations

Indicate the fortran include path and libraries for netcdf, e.g.:
```shell
export NC_INC="-I/usr/local/include"
export NC_LIB="-L/usr/local/lib -lnetcdf -lnetcdff"
```  

Use this script to create an ID for each ice shelf and attribute the corresponding melt rate from Rignot et al. (2013):
```shell
gfortran -c module_polarstereo.f90
gfortran -c $NC_INC identify_ice_shelf_mask_BEDMAP2.f90 
gfortran -o identify_ice_shelf_mask_BEDMAP2 module_polarstereo.o identify_ice_shelf_mask_BEDMAP2.o $NC_LIB
./identify_ice_shelf_mask_BEDMAP2
ls ice_shelf_mask_bedmap2_8km.nc  # created file
```

Use this script to interpolate the 3d observational thermal at the ice draft depths:
```shell
gfortran -c $NC_INC interpolate_TF_on_bedmap2_ice_draft.f90 
gfortran -o interpolate_TF_on_bedmap2_ice_draft  interpolate_TF_on_bedmap2_ice_draft.o $NC_LIB
./interpolate_TF_on_bedmap2_ice_draft
ls -al obs_thermal_forcing_on_bedmap2_ice_draft.nc  # created file
```

The following step is to calculate the parameters. The following options exist in calculate\_K0\_DeltaT\_quadratic.f90, the options used for ISMIP6 are in bold :
* para :
  + **= 1  quadratic non-local :  TF.<TF>**
  + **= 2  quadratic local : TF.TF**
* nn\_data\_melt :
  + = 1 using only Rignot et al. (2013)
  + = 2 using only Depoorter et al. (2013)
  + **= 3 using both Rignot et al. (2013) and Depoorter et al. (2013)**
* Nstat (number of samples to samples) :
  + = 10000 for quick tests
  + **= 100000 for more reproducible digits**
* Nsmooth (size of the running window, running mean between -Nsmooth/2 and +Nsmooth/2) : 
  + = 5
  + **= 0**
* ll\_upscaling :
  + = .true. to consider that BEDMAP2 represents the small glaciers that are part of the "upscaling" in Rignot and Depoorter's papers. 
  + **= .false. to consider that our resolution only represents the monitored glaciers (>~100km2) [false recommended for 8km resolution].**
* ll\_TS\_uncertainty :
  + = .false. for not considering uncertainties on observational TF.
  + **= .true. to take into account uncertainties on observational TF.**
```shell
  gfortran -c $NC_INC calculate_K0_DeltaT_quadratic.f90
  gfortran -o calculate_K0_DeltaT_quadratic  calculate_K0_DeltaT_quadratic.o $NC_LIB
  ./calculate_K0_DeltaT_quadratic
  ls -lrt coeff_K0_DeltaT_*.nc  # created files 
```

Then, we readjust the deltaT values to get the correct present-day melt rate in each basin with the 5th, 50th and 95th of gamma0. For the non-local parameterization, the correction is relatively small (but not zero because of the asymmetry induced by accounting for uncertainties on observational TF), but for the local parameterization, it is the core of the calculation:

For the non-local parameterization:
```shell
  gfortran -c $NC_INC readjust_deltaT_non_local_and_save_melt_MEDIAN.f90
  gfortran -c $NC_INC readjust_deltaT_non_local_and_save_melt_5TH_PCT.f90 
  gfortran -c $NC_INC readjust_deltaT_non_local_and_save_melt_95TH_PCT.f90 
  gfortran -o readjust_deltaT_non_local_and_save_melt_MEDIAN  readjust_deltaT_non_local_and_save_melt_MEDIAN.o $NC_LIB
  gfortran -o readjust_deltaT_non_local_and_save_melt_5TH_PCT  readjust_deltaT_non_local_and_save_melt_5TH_PCT.o $NC_LIB
  gfortran -o readjust_deltaT_non_local_and_save_melt_95TH_PCT  readjust_deltaT_non_local_and_save_melt_95TH_PCT.o $NC_LIB
  ./readjust_deltaT_non_local_and_save_melt_MEDIAN
  ./readjust_deltaT_non_local_and_save_melt_5TH_PCT
  ./readjust_deltaT_non_local_and_save_melt_95TH_PCT
```

For the local parameterization:
```shell
  gfortran -c $NC_INC readjust_deltaT_local_and_save_melt_MEDIAN.f90
  gfortran -c $NC_INC readjust_deltaT_local_and_save_melt_5TH_PCT.f90
  gfortran -c $NC_INC readjust_deltaT_local_and_save_melt_95TH_PCT.f90
  gfortran -o readjust_deltaT_local_and_save_melt_MEDIAN  readjust_deltaT_local_and_save_melt_MEDIAN.o $NC_LIB
  gfortran -o readjust_deltaT_local_and_save_melt_5TH_PCT  readjust_deltaT_local_and_save_melt_5TH_PCT.o $NC_LIB
  gfortran -o readjust_deltaT_local_and_save_melt_95TH_PCT  readjust_deltaT_local_and_save_melt_95TH_PCT.o $NC_LIB
  ./readjust_deltaT_local_and_save_melt_MEDIAN
  ./readjust_deltaT_local_and_save_melt_5TH_PCT
  ./readjust_deltaT_local_and_save_melt_95TH_PCT
```

To finalize the files for ISMIP6:

```shell
gfortran -c $NC_INC rewrite_adjusted_deltaT_NON_LOCAL.f90
gfortran -o rewrite_adjusted_deltaT_NON_LOCAL rewrite_adjusted_deltaT_NON_LOCAL.o $NC_LIB
./rewrite_adjusted_deltaT_NON_LOCAL
./finalize_files_for_ISMIP6_NON_LOCAL.sh
```

```shell
  gfortran -c $NC_INC rewrite_adjusted_deltaT_LOCAL.f90
  gfortran -o rewrite_adjusted_deltaT_LOCAL rewrite_adjusted_deltaT_LOCAL.o $NC_LIB
  ./rewrite_adjusted_deltaT_LOCAL
  ./finalize_files_for_ISMIP6_LOCAL.sh
```
