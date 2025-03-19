#!/bin/bash

#SBATCH --job-name=EPIC-RT
  # %x gives job-name (SLURM_JOB_NAME)
  # %j gives jobid (individual SLURM_JOB_ID)
  # %A gives jobid (master     SLURM_ARRAY_JOB_ID)
  # %a gives array task id number
  #  https://slurm.schedmd.com/sbatch.html
#SBATCH --nodes=NODES
#SBATCH --ntasks=NTASKS
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --account=e710
#SBATCH --partition=standard
#SBATCH --qos=standard

# refer to ARCHER2 documentation for more information
#     https://docs.archer2.ac.uk/user-guide/

# set number of threads per process
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export FI_OFI_RXM_SAR_LIMIT=64K
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
export MPI_DIR=$MPICH_DIR
export NETCDF_C_DIR=$NETCDF_DIR
export NETCDF_FORTRAN_DIR=$NETCDF_DIR
export FC=ftn

# assumes MPI-enabled build

# adjust directory to top-level epic install directory (prefix)
EPIC_DIR="/work/e710/e710/mf248/gcc/epic"

# adjust config, assumes the case was prepared (.nc file present and namelists/config edited)
srun --unbuffered --distribution=block:block --hint=nomultithread ${EPIC_DIR}/bin/epic3d --config rt_epic_late_restart.config --restart epic_rt_NXxNYxNZ_prepare_0000000003_parcels.nc --verbose
