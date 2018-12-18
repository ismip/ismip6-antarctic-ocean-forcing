#!/bin/bash -l
# comment out if using debug queue
#SBATCH --partition=regular
# comment in to get premium queue
##SBATCH --qos=premium
# comment in to get the debug queue
##SBATCH --partition=debug
# change number of nodes to change the number of parallel tasks
# (anything between 1 and the total number of tasks to run)
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --account=m1041
#SBATCH --job-name=ismip6_@year
#SBATCH --output=ismip6_@year.o%j
#SBATCH --error=ismip6_@year.e%j
#SBATCH -L SCRATCH
#SBATCH -C haswell

export OMP_NUM_THREADS=1

module unload python python/base e3sm-unified
source /global/homes/x/xylar/miniconda2/etc/profile.d/conda.sh
conda activate ismip6_ocean_forcing
export HDF5_USE_FILE_LOCKING=FALSE

srun -N 1 -n 1 -c 1 python -m ismip6_ocean_forcing ccsm4/@year/config.ccsm4
