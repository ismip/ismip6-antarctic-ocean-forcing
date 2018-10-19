# ismip6-ocean-forcing

Scripts for generating ocean forcing (primarily Antarctic) for the ISMIP6 activity

## Required conda environment

To generate the ocean forcing data, a conda environment can be set up with
the required packages as follows:

```bash
conda create -n ismip6_ocean_forcing -c conda-forge python=3.6 numpy scipy \
    matplotlib netCDF4 xarray progressbar2 basemap descartes cartopy shapely \
    nco gsw scikit-fmm phshp
```

If you need help installing `conda`,
[see instructions here](https://docs.anaconda.com/anaconda/install/)

To use this environment, you simply run:
```bash
conda activate ismip6_ocean_forcing
```
Note: if the `conda` command is not found, conda was not added to your
`.bashrc` and you may need to read the `conda` documentation to figure out how
to activate the base conda environment.

