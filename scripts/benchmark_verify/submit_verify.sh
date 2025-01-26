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

export WORK_DIR=/work/e710/e710/mf248
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
    module load cray-dsmml/0.2.2
    module load cray-openshmemx/11.5.7
    module load cray-netcdf-hdf5parallel/4.9.0.1

    module load cpe/23.09

    export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_DIR/crayclang/14.0
    export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_DIR/crayclang/14.0
    export FC=ftn
    export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
fi

if test "GRAPH_TYPE" = "shmem"; then
    echo "Setting SHMEM symmetric size"
    export SHMEM_SYMMETRIC_SIZE=1G
    export SHMEM_VERSION_DISPLAY=0
    export SHMEM_ENV_DISPLAY=0
fi

export SLURM_CPU_FREQ_REQ=2000000

#conda activate epic-env
source /work/e710/e710/mf248/miniconda3/bin/activate epic-env
export EXEC_PATH=/work/e710/e710/mf248/COMPILER/clustering/bin/pytools
export PYTHONPATH=$PYTHONPATH:$WORK_DIR/COMPILER/clustering/bin/pytools

PATH=/work/e710/e710/mf248/COMPILER/clustering/bin:$PATH

if test "GRAPH_TYPE" = "shmem"; then
    echo "Run OpenSHMEM"
    python ${EXEC_PATH}/verify_cluster_algorithm.py \
	    --n_ranks 16 32 64 128 256 \
	    --n_parcel_per_cell 40 \
	    --nx 32 \
	    --ny 32 \
	    --nz 32 \
	    --min_vratio 40.0 \
	    --verbose \
 	    --n_samples N_SAMPLES \
	    --cmd srun \
	    --seed SEED \
    	    --graph-type "GRAPH_TYPE"
else
    echo "Run GRAPH_TYPE"
    python ${EXEC_PATH}/verify_cluster_algorithm.py \
            --n_ranks 16 32 64 128 256 \
            --n_parcel_per_cell 40 \
            --nx 32 \
            --ny 32 \
            --nz 32 \
            --min_vratio 40.0 \
            --verbose \
            --n_samples N_SAMPLES \
            --cmd srun \
            --seed SEED \
            --graph-type "GRAPH_TYPE" \
	    --subcomm
fi
