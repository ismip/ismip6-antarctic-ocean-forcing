Download and Process MIROC-ESM-CHEM thetao and so
=================================================

A CERA or PCMDI account account is required to download the MIROC-ESM-CHEM CMIP5
data sets from which the forcing time series is derived.  Log in to your
account, then download the zip file ("Download all files") from the following
four links:

[MIROC-ESM-CHEM historical run so](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=MIM7hiMOOso111v111004)

[MIROC-ESM-CHEM historical run thetao](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=MIM7hiMOOthetao111v111004)

[MIROC-ESM-CHEM RCP8.5 run so](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=MIM7r8MOOso111v111004)

[MIROC-ESM-CHEM RCP8.5 run thetao](https://cera-www.dkrz.de/WDCC/ui/cerasearch/downloadForm?acronym=MIM7r8MOOthetao111v111004)

Finally, run `process_miroc_esm_chem.py -o /path/to/miroc_esm_chem/files/` with
the path where the MIROC-ESM-CHEM data sets are stored.

