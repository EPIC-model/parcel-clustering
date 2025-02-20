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

machine=''

# number of samples
n_samples=10

# RNG seed
seed=42
while getopts "h?m:n:s:": option; do
    case "$option" in
        h|\?)
            print_help
            exit 0
            ;;
        m)
            machine=$OPTARG
            ;;
        n)
            n_samples=$OPTARG
            ;;
        s)
            seed=$OPTARG
            ;;
    esac
done

if ! test "$machine" = "archer2" && ! test "$machine" = "cirrus"; then
    echo "Only 'archer2' and 'cirrus' machines supported. Exiting."
    exit 1
fi


# bin directories of executables:
source "../$machine.sh"
# --------------------------------------------------------


if ! test "$CONDA_EXE"; then
    echo "No CONDA environment. Checking if 'python_exe' is set."

    if ! test "$python_exe"; then
        exit 1
    fi

    CONDA_EXE=$python_exe
fi

j=0
for bin_dir in ${bins[*]}; do
    compiler="${compilers[$j]}"
    with_caf="${enable_caf[$j]}"

    # Coarray Fortran (CAF) is a separate build:
    if test "$with_caf" = "yes"; then
        run_job $machine $compiler "$bin_dir" "caf" $n_samples $seed $conda_env
    else
        for i in "p2p" "rma" "shmem"; do
            run_job $machine $compiler "$bin_dir" $i $n_samples $seed $conda_env
        done
    fi

    j=$((j+1))
done
