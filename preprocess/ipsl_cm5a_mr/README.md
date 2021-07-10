Download and Process IPSL-CM5A-MR thetao and so
===============================================

An [ESGF account](https://esgf-node.llnl.gov/login/) account is required to
download the IPSL-CM5A-MR CMIP5 data sets from which the forcing time series is
derived.  Log in to your account, then search for model: IPSL-CM5A-MR,
ensemble: r1i1p1, expereiment: historical,rcp85, variable: so,thetao.
Turn on "Show Replicas" in case the main site is down (happens relatively
often).

Either choose "List Files" or "WGET Script".  If the former, find files with
so and thetao data.  If the latter, edit the script to only download the files
needed.  In either case, from the historical, only files from 1994 on are
needed.

Finally, run:
```bash
./process_ipsl_cm5a_mr_rcp85.py -o /path/to/ipsl_cm5a_mr/files/
```
with the path where the IPSL-CM5A-MR RCP8.5 data sets are stored and similarly
for RCP2.6.

