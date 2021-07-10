# Instructions for processing obs and model data into ISMIP6 format

I process the CCSM4 data on one of the supercomputers I have access to.  I run
each year as a separate job (each takes a couple of hours running on 8 threads
on a single node -- I run out of memory with more threads than 8).  You could
process many years in a single job if you prefer, it would just take longer.

Then, once all years are processed, I combine groups of 20 years into a single
file, then I combine all of those files into a single long time series.

Here would be the approximate procedure, based on the scripts and templates in
this directory:

1. Clone the https://github.com/ismip/ismip6-antarctic-ocean-forcing repo
   somewhere

2. Setup the conda environment you need:
```bash
conda create -n ismip6_ocean_forcing -c conda-forge --file dev-spec.txt

python -m pip install -e .

conda activate ismip6_ocean_forcing
```

3. Download or copy `config.obs` and `job_script_obs.bash` from this directory
   into the same directory where you put the symlink

4. Modify config.obs as needed (mainly to change the number of parallel tasks
   -- threads -- you can afford to run on a single node of your machine)

5. Modify and submit `job_script_obs.bash` or just run the python command if
   you're running on a machine where you don't need to submit jobs:
```bash
python -m ismip6_ocean_forcing config.obs
```

7. Make a directory for your GCM data (e.g. `cosmos`)

8. Download/copy all the other files in this directory into the new directory

9. If you are working with a model already suppored by `ismip6_ocean_forcing`,
   go into the `preprocess` directory in this repo, find the model, and follow
   the steps outlined there to download and preprocess the data.

10. If you are using a different model, you will need to preprocess the data
    in a similar fashion, resulting in similar file names.

11. Modify the templates to suit your needs (e.g. different job scripts for
    whatever machine you use).   For exmaple, replace `ccsm4` with `cosmos` and
    `CCSM4` with `COSMOS`. In config.template, the options you would most likely
    change would be:

    ```
    # The name of the model
    name = CCSM4

    # Input files
    temperatureFileName = ccsm4/thetao_annual_CCSM4_rcp85_r1i1p1_185001-210012.nc
    salinityFileName = ccsm4/so_annual_CCSM4_rcp85_r1i1p1_185001-210012.nc
    
    # Variable names
    lon = lon
    lat = lat
    z = lev
    time = time
    temperature = thetao
    salinity = so
    z_bnds = lev_bnds
    ```

12. Run `setup_years.py` to make a folder for each year (or modify the script
    to process multiple years at a time)

13. Submit just one job for one year (e.g. `1995`) so some common data gets
    processed (mapping files to the 8km grid and matrices for doing
    extrapolation/interpolation under ice shelves).

14. Once this job finishes, you can submit all the jobs for all the other
    years.  (They share mapping files and matrices, so that is the reason to
    not run them all at once.)

15. Run `setup_decades.py` to make a set of folders with job scripts and config
    files needed for combining and processing blocks of 20 years of data

16. Once all years have been processed individually, submit one job for one of
    the 20-year periods (doesn't matter which).  This will compute a
    climatology from which anomalies are computed, then a time series for the
    20-year block of time you chose, then anomalies with respect to the
    climatology, then will add the observational climatology to the anomalies.

17. Once the first 20-year data set has been processed, you can submit jobs for
    all the remaining 20-year data sets.  (They use the same reference
    climatology, so that is the reason to not run them all at once.)

18. Once all the 20-year data sets have been processed, you can combine them
    into one (or two -- historical and rcp8.5 in my example) final time series
    by modifying and running `combine_all.py`

Obviously, you will probably need to make tweaks to these files along the way.
You may also find that you run into bugs in my python code, in which case I
would very much like you to submit and issue:
https://github.com/ismip/ismip6-antarctic-ocean-forcing/issues
or a pull request with a bug fix.
