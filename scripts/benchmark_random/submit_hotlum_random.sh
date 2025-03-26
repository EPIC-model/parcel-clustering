#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --output=%x.o%j
#SBATCH --time=01:00:00
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
# Set the number of threads to 1
#   This prevents any threaded system libraries from automatically
#   using threading.
export OMP_NUM_THREADS=1
export OMP_PLACES=cores

# Set the eager limit
# (see also https://docs.archer2.ac.uk/user-guide/tuning/#setting-the-eager-limit-on-archer2)
export FI_OFI_RXM_SAR_LIMIT=64K

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

export MPLCONFIGDIR=$PWD

if test "COMPILER" = "gnu"; then
    echo "Loading the GNU Compiler Collection (GCC)"
    module switch PrgEnv-cray PrgEnv-gnu
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load cray-dsmml
    module load cray-openshmemx
    module load cpe/24.11
    export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
    export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
    export FC=ftn
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
elif test "COMPILER" = "cray"; then
    echo "Loading the Cray Compiling Environment (CCE)"
    module load PrgEnv-cray
    module load cce
    module load cray-mpich
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load cray-dsmml
    module load cray-openshmemx

    # load latest modules
    module load cpe/24.11
    export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
    export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
    export FC=ftn
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
fi

#export MPIR_CVAR_CH4_OFI_ENABLE_RMA=0
echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS tasks."

module list

bin_dir=BIN_DIR
PATH=${bin_dir}:$PATH

rm /tmp/benchmark_random

sbcast --compress=none ${bin_dir}/benchmark_random /tmp/benchmark_random
for i in $(seq 1 NREPEAT); do
    srun -v --nodes=NODES \
        --ntasks=NTASKS \
        --unbuffered \
        --distribution=block:block \
	--hint=nomultithread \
        /tmp/benchmark_random \
        --nx NX \
        --ny NY \
        --nz NZ \
        --lx LX \
        --ly LY \
        --lz LZ \
        --min-vratio 20.0 \
        --nppc 20 \
        --niter NITER \
        --small-parcel-fraction SMALL_PARCEL_FRACTION \
        --shuffle \
        --csvfname "MACHINE-COMPILER-shmem-random-nx-NX-ny-NY-nz-NZ-nodes-NODES" \
        --comm-type "shmem"
    for g in "p2p" "rma"; do
        srun --nodes=NODES \
            --ntasks=NTASKS \
            --unbuffered \
            --distribution=block:block \
            --hint=nomultithread \
            /tmp/benchmark_random \
            --nx NX \
            --ny NY \
            --nz NZ \
            --lx LX \
            --ly LY \
            --lz LZ \
            --min-vratio 20.0 \
            --nppc 20 \
            --niter NITER \
            --small-parcel-fraction SMALL_PARCEL_FRACTION \
            --shuffle \
            --csvfname "MACHINE-COMPILER-$g-random-nx-NX-ny-NY-nz-NZ-nodes-NODES" \
            --comm-type "$g"

        if test "SUBCOMM" = "true"; then
            srun --nodes=NODES \
                --ntasks=NTASKS \
                --unbuffered \
                --distribution=block:block \
                --hint=nomultithread \
                /tmp/benchmark_random \
                --nx NX \
                --ny NY \
                --nz NZ \
                --lx LX \
                --ly LY \
                --lz LZ \
                --min-vratio 20.0 \
                --nppc 20 \
                --niter NITER \
                --small-parcel-fraction SMALL_PARCEL_FRACTION \
                --shuffle \
                --csvfname "MACHINE-COMPILER-$g-random-nx-NX-ny-NY-nz-NZ-nodes-NODES-subcomm" \
                --comm-type "$g" \
                --subcomm
        fi
    done
done

