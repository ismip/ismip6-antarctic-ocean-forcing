#!/bin/bash

rm -f coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc

ncks -F -d pct,1 coeff_K0_DeltaT_quadratic_local_Rignot_Depoorter_with_TS_uncert_ADJUSTED.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncwa -O -F -a pct,1 coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncks -O -x -v DeltaT_basin_pc05,DeltaT_basin_pc95,pct coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncrename -O -v DeltaT_basin_pc50,deltaT_basin -v K0,gamma0 coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncatted -O -a units,deltaT_basin,m,c,'K' coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncatted -O -a units,gamma0,m,c,'m/yr' coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncatted -O -a long_name,gamma0,m,c,'gamma_0 coefficient' coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncatted -O -a cell_methods,,d,, coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc
ncatted -h -O -a ,global,d,, coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_5th_percentile.nc

ncks -F -d pct,2 coeff_K0_DeltaT_quadratic_local_Rignot_Depoorter_with_TS_uncert_ADJUSTED.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncwa -O -F -a pct,1 coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncks -O -x -v DeltaT_basin_pc05,DeltaT_basin_pc95,pct coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncrename -O -v DeltaT_basin_pc50,deltaT_basin -v K0,gamma0 coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncatted -O -a units,deltaT_basin,m,c,'K' coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncatted -O -a units,gamma0,m,c,'m/yr' coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncatted -O -a long_name,gamma0,m,c,'gamma_0 coefficient' coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncatted -O -a cell_methods,,d,, coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc
ncatted -h -O -a ,global,d,, coeff_gamma0_DeltaT_quadratic_local_median.nc coeff_gamma0_DeltaT_quadratic_local_median.nc

ncks -F -d pct,3 coeff_K0_DeltaT_quadratic_local_Rignot_Depoorter_with_TS_uncert_ADJUSTED.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncwa -O -F -a pct,1 coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncks -O -x -v DeltaT_basin_pc05,DeltaT_basin_pc95,pct coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncrename -O -v DeltaT_basin_pc50,deltaT_basin -v K0,gamma0 coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncatted -O -a units,deltaT_basin,m,c,'K' coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncatted -O -a units,gamma0,m,c,'m/yr' coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncatted -O -a long_name,gamma0,m,c,'gamma_0 coefficient' coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncatted -O -a cell_methods,,d,, coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
ncatted -h -O -a ,global,d,, coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc coeff_gamma0_DeltaT_quadratic_local_95th_percentile.nc
