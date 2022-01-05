Download and Process UKESM1-0-LL thetao and so
==============================================

An Earth System Grid account is required to download the UKESM1-0-LL
CMIP6 data sets from which the forcing time series is derived.  Log
in to your account, go to CMIP6, select UKESM1-0-LL as the
Source ID; historical and ssp126 as the Experiment IC; r4i1p1f2 as
the Variant Label; and so and thetao in Variables.  Click "Search".  You many
want to choose "Show All Replicas" and "Search" again.

Download the WGET Script for each Source IC and Variable, and run the scripts.
If all goes correctly, you should only need to log in once with your OpenID to
run all the WGET Scripts.

Finally, run `process_ukesm1_0_ll.py -o /path/to/ukesm1-0-ll/files/` with the
path where the UKESM1-0-LL data sets are stored.

