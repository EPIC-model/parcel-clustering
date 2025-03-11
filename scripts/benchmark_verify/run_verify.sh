#!/bin/bash

run_job() {

    local machine=${1}
    local fname="submit_${machine}_verify.sh"

    local compiler=${2}
    local bin_dir=${3}
    local comm_type=${4}
    local n_samples=${5}
    local seed=${6}

    echo "--------------------------------"
    echo "Run jobs with following options:"
    echo "machine    = $machine"
    echo "compiler   = $compiler"
    echo "bin_dir    = $bin_dir"
    echo "comm_type  = $comm_type"
    echo "n_samples  = $n_samples"
    echo "seed       = $seed"
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

    echo "Submit job $comm_type with $compiler. Running $n_samples samples with seed $seed."
    sbatch $fname
    cd "../.."
}

print_help() {
    echo "This script submits code verification jobs."
    echo "Arguments:"
    echo "    -h    print this message"
    echo "    -c    communication layer, valid options: 'p2p', 'shmem' and 'rma'"
    echo "    -m    machine to run on, e.g. 'archer2', 'cirrus', 'hotlum'"
    echo "          (requirement: <machine>.sh and 'submit_<machine>_verify.sh)"
    echo "    -n    number of random samples"
    echo "    -s    seed for RNG"
}

check_for_input() {
    if ! test "${2}"; then
        echo "Please specify '${1}'. Exiting."
	exit 1
    fi
}

# --------------------------------------------------------
# User options:

machine=''

# number of samples
n_samples=10

# RNG seed
seed=42

comm="p2p"

while getopts "h?m:n:s:c:": option; do
    case "$option" in
	c)
	    comm=$OPTARG
	    ;;
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

check_for_input "machine" $machine
check_for_input "number of samples" $n_samples
check_for_input "seed" $seed
check_for_input "communication layer" $comm

l_comm_valid=0
for i in "p2p" "rma" "shmem"; do
    if test "$comm" = "$i"; then
        l_comm_valid=1
    fi
done

if ! test $l_comm_valid = 1; then
    echo "Invalid communication layer. Exiting."
    exit 1
fi

if ! test -f "../$machine.sh"; then
    echo "Unable to run on $machine. The file ${machine}.sh does not exist. Exiting."
    exit 1
fi

if ! test -f "submit_${machine}_verify.sh" ; then
    echo "Unable to run on $machine. The file submit_${machine}_verify.sh does not exist. Exiting."
    exit 1
fi


# bin directories of executables:
source "../$machine.sh"
# --------------------------------------------------------

j=0
for bin_dir in ${bins[*]}; do
    compiler="${compilers[$j]}"

    run_job $machine $compiler "$bin_dir" $comm $n_samples $seed

    j=$((j+1))
done
