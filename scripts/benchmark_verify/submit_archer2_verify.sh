#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --output=%x.o%j
#SBATCH --time=96:00:00
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --constraint=StandardMem
#SBATCH --switches=1
#SBATCH --account=e710
#SBATCH --partition=standard
#SBATCH --qos=long

# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export FI_OFI_RXM_SAR_LIMIT=64K
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

if test "COMPILER" = "gnu"; then
    echo "Loading the GNU Compiler Collection (GCC)"
    module load PrgEnv-gnu
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load cray-dsmml
    module load cray-openshmemx
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_NETCDF_HDF5PARALLEL_PREFIX/lib
elif test "COMPILER" = "cray"; then
    echo "Loading the Cray Compiling Environment (CCE)"
    module load PrgEnv-cray/8.3.3
    module load cce/15.0.0
    module load cray-mpich/8.1.23
    module load cray-hdf5-parallel/1.12.2.1
    module load cray-dsmml/0.2.2
    module load cray-openshmemx/11.5.7
    module load cray-netcdf-hdf5parallel/4.9.0.1

    module load cpe/23.09

    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_NETCDF_HDF5PARALLEL_PREFIX/lib
fi

if test "COMM_TYPE" = "shmem"; then
    echo "Setting SHMEM symmetric size"
    export SHMEM_SYMMETRIC_SIZE=1G
    export SHMEM_VERSION_DISPLAY=0
    export SHMEM_ENV_DISPLAY=0
fi

export SLURM_CPU_FREQ_REQ=2000000

module list

bin_dir=BIN_DIR

PATH=${bin_dir}:$PATH

echo "Run COMM_TYPE"

if test "COMM_TYPE" = "shmem"; then
    ${bin_dir}/verify_cluster_algorithm \
        --nranks 32 64 96 128 160 192 224 256 \
        --ntasks-per-node 128 \
        --nppc 40 \
        --nx 32 \
        --ny 32 \
        --nz 32 \
        --min-vratio 40.0 \
        --verbose \
        --nsamples N_SAMPLES \
        --cmd srun \
        --seed SEED \
        --comm-type "COMM_TYPE"
else
    ${bin_dir}/verify_cluster_algorithm \
        --nranks 32 64 96 128 160 192 224 256 \
        --ntasks-per-node 128 \
        --nppc 40 \
        --nx 32 \
        --ny 32 \
        --nz 32 \
        --min-vratio 40.0 \
        --verbose \
        --nsamples N_SAMPLES \
        --cmd srun \
        --seed SEED \
        --comm-type "COMM_TYPE" \
        --subcomm
fi
