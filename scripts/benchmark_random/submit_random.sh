#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --output=%x.o%j
#SBATCH --time=00:10:00
#SBATCH --nodes=NODES
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1
#SBATCH --constraint=StandardMem
#SBATCH --switches=1
#SBATCH --account=e710
#SBATCH --partition=standard
#SBATCH --qos=standard

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

    # make gcc/12.2.0 available and load it
    module load load-epcc-module;
    module load  extra-compilers/1.0

    # update all other modules:
    module load cpe/23.09

    export NETCDF_C_DIR=$NETCDF_DIR
    export NETCDF_FORTRAN_DIR=$NETCDF_DIR
    export FC=ftn
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
elif test "COMPILER" = "cray"; then
    echo "Loading the Cray Compiling Environment (CCE)"
    module load PrgEnv-cray/8.3.3
    module load cce/15.0.0
    module load cray-mpich/8.1.23
    module load cray-hdf5-parallel/1.12.2.1
    module load cray-netcdf-hdf5parallel/4.9.0.1
    module load cray-dsmml/0.2.2
    module load cray-openshmemx/11.5.7

    # load latest modules
    module load cpe/23.09
    export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_DIR/crayclang/14.0
    export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_DIR/crayclang/14.0
    export FC=ftn
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
fi

#export PGAS_MEMINFO_DISPLAY=1
#export XT_SYMMETRIC_HEAP_SIZE=1G
#export CRAY_PGAS_MAX_CONCURRENCY=1

echo "Setting SHMEM symmetric size"
export SHMEM_SYMMETRIC_SIZE=1G
export SHMEM_VERSION_DISPLAY=0
export SHMEM_ENV_DISPLAY=0

export SLURM_CPU_FREQ_REQ=2000000

module load craype-hugepages2M
export HUGETLB_VERBOSE=2

echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS tasks."

bin_dir=BIN_DIR
PATH=${bin_dir}:$PATH

sbcast --compress=none ${bin_dir}/benchmark_random /tmp/benchmark_random
for i in $(seq 1 NREPEAT); do
    if test "ENABLE_CAF" = "yes"; then
        srun --nodes=NODES \
            --ntasks=NTASKS \
            --unbuffered \
            --distribution=block:block \
            /tmp/benchmark_random \
            --nx NX \
            --ny NY \
            --nz 32 \
            --lx LX \
            --ly LY \
            --lz 10.0 \
            --xlen LX \
            --ylen LY \
            --zlen 10.0 \
            --min_vratio 20.0 \
            --n_per_cell 20 \
            --niter NITER \
            --shuffle \
            --ncfname "COMPILER-caf-random-nx-NX-ny-NY-nodes-NODES.nc" \
            --graph-type "caf"
    else
        srun --nodes=NODES \
            --ntasks=NTASKS \
            --unbuffered \
            --distribution=block:block \
            /tmp/benchmark_random \
            --nx NX \
            --ny NY \
            --nz 32 \
            --lx LX \
            --ly LY \
            --lz 10.0 \
            --xlen LX \
            --ylen LY \
            --zlen 10.0 \
            --min_vratio 20.0 \
            --n_per_cell 20 \
            --niter NITER \
            --shuffle \
            --ncfname "COMPILER-shmem-random-nx-NX-ny-NY-nodes-NODES.nc" \
            --graph-type "shmem"
        for g in "p2p" "rma"; do
            srun --nodes=NODES \
                --ntasks=NTASKS \
                --unbuffered \
                --distribution=block:block \
                --hint=nomultithread \
                /tmp/benchmark_random \
                --nx NX \
                --ny NY \
                --nz 32 \
                --lx LX \
                --ly LY \
                --lz 10.0 \
                --xlen LX \
                --ylen LY \
                --zlen 10.0 \
                --min_vratio 20.0 \
                --n_per_cell 20 \
                --niter NITER \
                --shuffle \
                --ncfname "COMPILER-$g-random-nx-NX-ny-NY-nodes-NODES.nc" \
                --graph-type "$g"

            if test "SUBCOMM" = "true"; then
                srun --nodes=NODES \
                    --ntasks=NTASKS \
                    --unbuffered \
                    --distribution=block:block \
                    --hint=nomultithread \
                    /tmp/benchmark_random \
                    --nx NX \
                    --ny NY \
                    --nz 32 \
                    --lx LX \
                    --ly LY \
                    --lz 10.0 \
                    --xlen LX \
                    --ylen LY \
                    --zlen 10.0 \
                    --min_vratio 20.0 \
                    --n_per_cell 20 \
                    --niter NITER \
                    --shuffle \
                    --ncfname "COMPILER-$g-random-nx-NX-ny-NY-nodes-NODES-subcomm.nc" \
                    --graph-type "$g" \
                    --subcomm
            fi
        done
    fi
done

