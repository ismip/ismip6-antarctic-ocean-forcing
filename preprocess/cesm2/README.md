Download and Process CESM2 thetao and so
========================================

An Earth System Grid account is required to download the CESM2
CMIP6 data sets from which the forcing time series is derived.  Log
in to your account, go to CMIP6, select CESM as the
Source ID; historical as the Experiment IC; r11i1p1f1 as
the Variant Label; and so and thetao in Variables.  Click "Search".  You many
want to choose "Show All Replicas" and "Search" again.

Download the WGET Script for the "gn" (native grid) verison of  each Variable.

Now, for the ssp data, go to CMIP6, select CESM as the
Source ID; ssp585 as the Experiment IC; r2i1p1f1 as
the Variant Label; and so and thetao in Variables.  Click "Search".  You many
want to choose "Show All Replicas" and "Search" again.

Again, download the WGET Script for the "gn" version of each Variable.  Run all
of the WGET scripts. If all goes correctly, you should only need to log in once
with your OpenID to run all the WGET Scripts.

Finally, run `process_cesm2_ssp585.py -o /path/to/cesm2/files/` with the
path where the CESM2 data sets are stored (or similarly for the "historical"
script).

