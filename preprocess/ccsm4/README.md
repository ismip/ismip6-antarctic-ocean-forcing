Download and Process CCSM4 thetao and so
========================================

An [ESGF account](https://esgf-node.llnl.gov/login/) account is required to
download the CCSM4 CMIP5 data sets from which the forcing time series is
derived.  Log in to your account, then search for model: `CCSM4`,
ensemble: `r1i1p1`, experiment: `historical,rcp85`, variable: `so,thetao`.
Turn on "Show Replicas" in case the main site is down (happens relatively
often).

Either choose "List Files" or "WGET Script".  If the former, find files with
`so` and `thetao` data.  If the latter, edit the script to only download the 
files  needed.

Finally, run:
```bash
./process_ccsm4.py -o /path/to/ccsm4/files/
```
with the path where the CCSM4 RCP8.5 data sets are stored.
