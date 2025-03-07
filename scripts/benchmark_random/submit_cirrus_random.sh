#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --output=%x.o%j
#SBATCH --time=00:10:00
#SBATCH --nodes=NODES
#SBATCH --tasks-per-node=36
#SBATCH --cpus-per-task=1
#SBATCH --switches=1
#SBATCH --account=e710
#SBATCH --partition=standard
#SBATCH --qos=standard # largescale
#SBATCH --exclusive
#SBATCH --distribution=block:block

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
    module load libtool/2.4.7
    module load gcc/10.2.0
    module load openmpi/4.1.6
    module load hdf5parallel/1.14.3-gcc10-ompi416
    export NETCDF_C_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
    export NETCDF_FORTRAN_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
elif test "COMPILER" = "intel-mpi"; then
    echo "Loading the Intel Compiler Environment"
    module load libtool/2.4.7
    module load intel-20.4/mpi
    module load intel-20.4/compilers
    module load netcdf-parallel/4.9.2-intel20-impi20
    export NETCDF_C_DIR=$NETCDF_DIR
    export NETCDF_FORTRAN_DIR=$NETCDF_DIR
elif test "COMPILER" = "intel-mpt"; then
    echo "Loading the HPE MPT Environment with the Intel compiler"
    module load mpt/2.25
    module load intel-20.4/compilers
    module load netcdf-parallel/4.9.2-intel20-mpt225
    export NETCDF_FORTRAN_DIR=$NETCDF_DIR
    export NETCDF_C_DIR=$NETCDF_DIR
    export MPICC_CC=icc
    export MPICXX_CXX=icpc
fi

echo "Setting SHMEM symmetric size"
export SHMEM_SYMMETRIC_SIZE=1G
export SHMEM_VERSION_DISPLAY=0
export SHMEM_ENV_DISPLAY=0

echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS tasks."

bin_dir=BIN_DIR
PATH=${bin_dir}:$PATH

for i in $(seq 1 NREPEAT); do
    if test "ENABLE_CAF" = "yes"; then
        srun --kill-on-bad-exit \
	     --nodes=NODES \
             --ntasks=NTASKS \
             --unbuffered \
             --distribution=block:block \
             ${bin_dir}/benchmark_random \
             --nx NX \
             --ny NY \
             --nz NZ \
             --lx LX \
             --ly LY \
             --lz LZ \
             --xlen LX \
             --ylen LY \
             --zlen LZ \
             --min-vratio 20.0 \
             --nppc 20 \
             --niter NITER \
             --shuffle \
             --ncfname "COMPILER-caf-random-nx-NX-ny-NY-nz-NZ-nodes-NODES.nc" \
             --comm-type "caf"
    else
        srun --kill-on-bad-exit \
	     --nodes=NODES \
             --ntasks=NTASKS \
             --unbuffered \
             --distribution=block:block \
             ${bin_dir}/benchmark_random \
             --nx NX \
             --ny NY \
             --nz NZ \
             --lx LX \
             --ly LY \
             --lz LZ \
             --xlen LX \
             --ylen LY \
             --zlen LZ \
             --min-vratio 20.0 \
             --nppc 20 \
             --niter NITER \
             --shuffle \
             --ncfname "COMPILER-shmem-random-nx-NX-ny-NY-nz-NZ-nodes-NODES.nc" \
             --comm-type "shmem"
        for g in "p2p" "rma"; do
            srun --kill-on-bad-exit \
	         --nodes=NODES \
                 --ntasks=NTASKS \
                 --unbuffered \
                 --distribution=block:block \
                 --hint=nomultithread \
                 ${bin_dir}/benchmark_random \
                 --nx NX \
                 --ny NY \
                 --nz NZ \
                 --lx LX \
                 --ly LY \
                 --lz LZ \
                 --xlen LX \
                 --ylen LY \
                 --zlen LZ \
                 --min-vratio 20.0 \
                 --nppc 20 \
                 --niter NITER \
                 --shuffle \
                 --ncfname "COMPILER-$g-random-nx-NX-ny-NY-nz-NZ-nodes-NODES.nc" \
                 --comm-type "$g"

            if test "SUBCOMM" = "true"; then
                srun --kill-on-bad-exit \
	             --nodes=NODES \
                     --ntasks=NTASKS \
                     --unbuffered \
                     --distribution=block:block \
                     --hint=nomultithread \
                     ${bin_dir}/benchmark_random \
                     --nx NX \
                     --ny NY \
                     --nz NZ \
                     --lx LX \
                     --ly LY \
                     --lz LZ \
                     --xlen LX \
                     --ylen LY \
                     --zlen LZ \
                     --min-vratio 20.0 \
                     --nppc 20 \
                     --niter NITER \
                     --shuffle \
                     --ncfname "COMPILER-$g-random-nx-NX-ny-NY-nz-NZ-nodes-NODES-subcomm.nc" \
                     --comm-type "$g" \
                     --subcomm
            fi
        done
    fi
done

