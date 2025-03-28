#!/bin/bash
#SBATCH --job-name=JOBNAME
#SBATCH --output=%x.o%j
#SBATCH --time=TIMELIMIT
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

if test "COMPILER" = "gnu"; then
    echo "Loading the GNU Compiler Collection (GCC)"
    module load libtool/2.4.7
    module load gcc/10.2.0
    export MPI_ROOT=/work/e710/e710/mf248/gcc/10.2.0/openmpi/5.0.7
    export PATH=$MPI_ROOT/bin:$PATH
    export NETCDF_C_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
    export NETCDF_FORTRAN_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
    export PATH=$PATH:$NETCDF_C_DIR/bin
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_C_DIR/lib:$MPI_ROOT/lib
else
    echo "Only GNU Compiler Collection available!"
    exit 1
fi

echo "Setting SHMEM symmetric size"
export SHMEM_SYMMETRIC_SIZE=1G
export SHMEM_VERSION_DISPLAY=0
export SHMEM_ENV_DISPLAY=0

echo "Running on $SLURM_NNODES nodes with $SLURM_NTASKS tasks."

module list

cmd=$(which mpirun)
echo "MPIRUN: $cmd"

bin_dir=BIN_DIR
PATH=${bin_dir}:$PATH

for i in $(seq 1 NREPEAT); do
    mpirun --bind-to core --mca spml ucx --mca osc ucx --mca pml ucx -np NTASKS \
         ${bin_dir}/benchmark_read \
         --dirname DIRNAME \
         --ncbasename NC_BASENAME \
         --niter NITER \
         --offset OFFSET \
         --nfiles NFILES \
         --size-factor SIZE_FACTOR \
         --csvfname "MACHINE-COMPILER-shmem-read-NAMETAG-nx-NX-ny-NY-nz-NZ-nodes-NODES" \
         --comm-type "shmem"
    for g in "p2p" "rma"; do
        mpirun --bind-to core --mca spml ucx --mca osc ucx --mca pml ucx -np NTASKS \
             ${bin_dir}/benchmark_read \
             --dirname DIRNAME \
             --ncbasename NC_BASENAME \
             --niter NITER \
             --offset OFFSET \
             --nfiles NFILES \
             --size-factor SIZE_FACTOR \
             --csvfname "MACHINE-COMPILER-$g-read-NAMETAG-nx-NX-ny-NY-nz-NZ-nodes-NODES" \
             --comm-type "$g"

        if test "SUBCOMM" = "true"; then
            mpirun --bind-to core --mca spml ucx --mca osc ucx --mca pml ucx -np NTASKS \
                 ${bin_dir}/benchmark_read \
                 --dirname DIRNAME \
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

