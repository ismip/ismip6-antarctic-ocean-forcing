Download and Process CCSM4 thetao and so
========================================

An Earth System Grid account is required to download the CCSM4 CMIP5 data sets
from which the forcing time series is derived.  Log in to your ESG account, then
download the following list of files for the historical data from the
[CCSM4 historical run](https://www.earthsystemgrid.org/dataset/cmip5.output1.NCAR.CCSM4.historical.mon.ocean.Omon.r1i1p1/file.html?filter=&variable=so&variable=thetao)

* `thetao_Omon_CCSM4_historical_r1i1p1_185001-185912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_186001-186912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_187001-187912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_188001-188912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_189001-189912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_190001-190912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_191001-191912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_192001-192912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_193001-193912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_194001-194912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_195001-195912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_196001-196912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_197001-197912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_198001-198912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_199001-199912.nc`
* `thetao_Omon_CCSM4_historical_r1i1p1_200001-200512.nc`
* `so_Omon_CCSM4_historical_r1i1p1_*.nc`

Next, download the following files from the
[CCSM4 RCP 8.5 run](https://www.earthsystemgrid.org/dataset/cmip5.output1.NCAR.CCSM4.rcp85.mon.ocean.Omon.r1i1p1/file.html?filter=&variable=so&variable=thetao)

* `thetao_Omon_CCSM4_rcp85_r1i1p1_200601-200912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_201001-201912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_202001-202912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_203001-203912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_204001-204912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_205001-205912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_206001-206912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_207001-207912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_208001-208912.nc`
* `thetao_Omon_CCSM4_rcp85_r1i1p1_209001-210012.nc`
* `so_Omon_CCSM4_rcp85_r1i1p1_*.nc`

Finally, run `process_ccsm4.py -o /path/to/ccsm4/files/` with the path where the
CCSM4 data sets are stored.

