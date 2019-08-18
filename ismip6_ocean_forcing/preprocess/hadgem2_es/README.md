Download and Process HadGEM2-ES thetao and so
=============================================

An [ESGF account](https://esgf-node.llnl.gov/login/) account is required to
download the HadGEM2-ES CMIP5 data sets from which the forcing time series is
derived.  Log in to your account, then search for model: HadGEM2-ES,
ensemble: r1i1p1, expereiment: historical,rcp85, variable: so,thetao.
Turn on "Show Replicas" in case the main site is down (happens relatively
often).

Either choose "List Files" or "WGET Script".  If the former, find files with
so and thetao data.  If the latter, edit the script to only download the files
needed.  In either case, from the historical, only files from 1994 on are
needed.

Finally, run:
```bash
./process_hadgem2_es_rcp85.py -o /path/to/hadgem2_es_rcp85/files/
```
with the path where the HadGEM2-ES RCP8.5 data sets are stored.

