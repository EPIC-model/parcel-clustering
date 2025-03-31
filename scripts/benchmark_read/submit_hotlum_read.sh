#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --output=%x.o%j
#SBATCH --time=00:10:00
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
    module load PrgEnv-gnu
    module load cray-hdf5-parallel
    module load cray-netcdf-hdf5parallel
    module load cray-dsmml
    module load cray-openshmemx
    module load cpe/24.11
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_NETCDF_HDF5PARALLEL_PREFIX/lib
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
    module load -f cpe/24.11
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CRAY_NETCDF_HDF5PARALLEL_PREFIX/lib
fi

echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS tasks."

module list

bin_dir=BIN_DIR
PATH=${bin_dir}:$PATH

sbcast --compress=none ${bin_dir}/benchmark_read /tmp/benchmark_read
for i in $(seq 1 NREPEAT); do
    srun --nodes=NODES \
        --ntasks=NTASKS \
        --unbuffered \
        --distribution=block:block \
        /tmp/benchmark_read \
        --ncbasename NC_BASENAME \
        --niter NITER \
        --offset OFFSET \
        --nfiles NFILES \
        --size-factor SIZE_FACTOR \
        --csvfname "MACHINE-COMPILER-$g-read-NAMETAG-nx-NX-ny-NY-nz-NZ-nodes-NODES" \
        --comm-type "shmem"
    for g in "p2p" "rma"; do
        srun --nodes=NODES \
            --ntasks=NTASKS \
            --unbuffered \
            --distribution=block:block \
            --hint=nomultithread \
            /tmp/benchmark_read \
            --ncbasename NC_BASENAME \
            --niter NITER \
            --offset OFFSET \
            --nfiles NFILES \
            --size-factor SIZE_FACTOR \
            --csvfname "MACHINE-COMPILER-$g-read-NAMETAG-nx-NX-ny-NY-nz-NZ-nodes-NODES" \
            --comm-type "$g"

        if test "SUBCOMM" = "true"; then
            srun --nodes=NODES \
                --ntasks=NTASKS \
                --unbuffered \
                --distribution=block:block \
                --hint=nomultithread \
                /tmp/benchmark_read \
                --ncbasename NC_BASENAME \
                --niter NITER \
                --offset OFFSET \
                --nfiles NFILES \
                --size-factor SIZE_FACTOR \
                --csvfname "MACHINE-COMPILER-$g-read-NAMETAG-nx-NX-ny-NY-nz-NZ-nodes-NODES-subcomm" \
                --comm-type "$g" \
                --subcomm
        fi
    done
done

