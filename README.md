# ismip6-antarctic-ocean-forcing

Scripts for generating Antarctic ocean forcing for the ISMIP6 activity

## Required conda environment

To generate the ocean forcing data, you need the `conda` package manager and 
the `mamba` dependency solver.  If you don't have a conda base environment with
`mamba` installed in it, please download and install
[Mambaforge](https://github.com/conda-forge/miniforge#mambaforge).

Then, this repository can be cloned and a conda  environment can be set up with
the required packages as follows:

```bash
git clone git@github.com:ismip/ismip6-antarctic-ocean-forcing.git
cd ismip6-antarctic-ocean-forcing
mamba create -n ismip6_ocean_forcing -c conda-forge --file dev-spec.txt
mamba activate ismip6_ocean_forcing
python -m pip install -e .
```

To use this environment, you simply run:
```bash
mamba activate ismip6_ocean_forcing
```
Note: if the `mamba` command is not found, `mamba` was not added to your
`.bashrc` or equivalent.  You will need to run something like:
```bash
source ~/mambaforge/etc/profile.d/conda.sh
source ~/mambaforge/etc/profile.d/mamba.sh
mamba activate
```
to ge the base environment with the `mamba command`.

