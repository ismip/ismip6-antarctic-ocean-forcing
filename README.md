# ismip6-antarctic-ocean-forcing

Scripts for generating Antarctic ocean forcing for the ISMIP6 activity

## Required conda environment

To generate the ocean forcing data, this repository can be cloned and a conda
environment can be set up with the required packages as follows:

```bash
git clone git@github.com:ismip/ismip6-antarctic-ocean-forcing.git
cd ismip6-antarctic-ocean-forcing
conda create -n ismip6_ocean_forcing -c conda-forge --file dev-spec.txt
conda activate ismip6_ocean_forcing
python -m pip install -e .
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

