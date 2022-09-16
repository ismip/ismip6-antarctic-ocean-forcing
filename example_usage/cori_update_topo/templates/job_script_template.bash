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
#SBATCH --account=m1795
#SBATCH --job-name=@modelLower_@year
#SBATCH --output=@modelLower/@year/ismip6.o%j
#SBATCH --error=@modelLower/@year/ismip6.e%j
#SBATCH -L SCRATCH
#SBATCH -C haswell

export OMP_NUM_THREADS=1

source ~/miniconda3/etc/profile.d/conda.sh
conda activate ismip6
export HDF5_USE_FILE_LOCKING=FALSE

srun -N 1 -n 1 -c 32 ismip6_ocean_forcing @modelLower/@year/config.@modelLower
