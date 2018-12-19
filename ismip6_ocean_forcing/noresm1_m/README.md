Download and Process NorESM1-M thetao and so
============================================

A CERA or PCMDI account account is required to download the NorESM1-M CMIP5
data sets from which the forcing time series is derived.  Log in to your
account, then download the zip file ("Download all files") from the following
four links:

[NorESM1-M historical run so](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=NCCNMhiMOOso111v110901)

[NorESM1-M historical run thetao](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=NCCNMhiMOOthetao111v110901)

[NorESM1-M RCP8.5 run so](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=NCCNMr8MOOso111v110901)

[NorESM1-M RCP8.5 run theta](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=NCCNMr8MOOthetao111v110901)

Finally, run `process_noresm1_m.py -o /path/to/noresm1_m/files/` with the path
where the NorESM1-M data sets are stored.

