#!/bin/bash

run_job() {

    local machine=${1}
    local fname="submit_${machine}_verify.sh"

    local compiler=${2}
    local bin_dir=${3}
    local comm_type=${4}
    local n_samples=${5}
    local seed=${6}
    local conda_dir=${CONDA_EXE%/*}
    local conda_env=${7}

    echo "--------------------------------"
    echo "Run jobs with following options:"
    echo "machine    = $machine"
    echo "compiler   = $compiler"
    echo "bin_dir    = $bin_dir"
    echo "comm_type  = $comm_type"
    echo "n_samples  = $n_samples"
    echo "seed       = $seed"
    echo "conda_dir  = $conda_dir"
    echo "conda_env  = $conda_env"
    echo "--------------------------------"

    mkdir -p "$compiler"
    cd "$compiler"

    mkdir -p "$comm_type"
    cd "$comm_type"
    cp "../../$fname" .

    sed -i "s:#SBATCH --job-name=JOBNAME:#SBATCH --job-name=$compiler-$comm_type:g" $fname
    sed -i "s:COMPILER:$compiler:g" $fname
    sed -i "s:COMM_TYPE:$comm_type:g" $fname
    sed -i "s:N_SAMPLES:$n_samples:g" $fname
    sed -i "s:SEED:$seed:g" $fname
    sed -i "s:BIN_DIR:$bin_dir:g" $fname
    sed -i "s:CONDA_DIR:$conda_dir:g" $fname
    sed -i "s:CONDA_ENV:$conda_env:g" $fname

    echo "Submit job $comm_type with $compiler. Running $n_samples samples with seed $seed."
    sbatch $fname
    cd "../.."
}

# --------------------------------------------------------
# User options:

# number of samples
n_samples=10

# RNG seed
seed=42

# PYTHON CONDA environment
conda_env="epic-env"

# bin directories of executables:
gnu_bin="/work/e710/e710/mf248/gnu/clustering/bin"
cray_bin="/work/e710/e710/mf248/cray/clustering/bin"
caf_bin="/work/e710/e710/mf248/cray-caf/clustering/bin"
# --------------------------------------------------------


if ! test "$CONDA_EXE"; then
    echo "No CONDA environment."
    exit 1
fi

for i in "p2p" "rma" "shmem"; do
    if test -d "$gnu_bin"; then
        run_job $machine "gnu" $gnu_bin $i $n_samples $seed $conda_env
    fi

    if test -d "$cray_bin"; then
        run_job $machine "cray" $cray_bin $i $n_samples $seed $conda_env
    fi
done

# Coarray Fortran (CAF) is a separate build:
if test -d "$caf_bin"; then
    run_job $machine "cray" $caf_bin "caf" $n_samples $seed $conda_env
fi
