Download and Process CNRM-CM6-1 and CNRM-ESM2-1 thetao and so
=============================================================

An Earth System Grid account is required to download the CNRM-CM6-1 and
CNRM-ESM2-1 CMIP6 data sets from which the forcing time series is derived.  Log
in to your account, go to CMIP6, select CNRM-CM6-1 or CNRM-ESM2-1 as the
Source ID; historical, ssp126 and/or ssp585 as the Experiment IC; r1i1p1f2 as
the Variant Label; and so and thetao in Variables.  Click "Search".  You many
want to choose "Show All Replicas" and "Search" again.

Download the WGET Script for each Source IC and Variable, and run the scripts.
If all goes correctly, you should only need to log in once with your OpenID to
run all the WGET Scripts.

Finally, run `process_cnrm_cm6_1.py -o /path/to/cnrm_cm6_1/files/` with the path
where the CNRM-CM6-1 data sets are stored, or similarly for CNRM-ESM-1

