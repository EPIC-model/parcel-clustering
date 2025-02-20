#!/bin/bash

run_jobs() {

    local machine=${1}
    local fname="submit_${machine}_random.sh"
    local ntasks_per_node=${2}

    local compiler=${3}

    local bin_dir=${4}
    local nrepeat=${5}
    local niter=${6}
    local nx=${7}
    local ny=${8}
    local nz=${9}
    local lx=${10}
    local ly=${11}
    local lz=${12}
    local min_ntasks=${13}
    local inc_ntasks=${14}
    local max_ntasks=${15}
    local subcomm=${16}
    local enable_caf=${17}

    echo "--------------------------------"
    echo "Run jobs with following options:"
    echo "machine         = $machine"
    echo "ntasks_per_node = $ntasks_per_node"
    echo "fname           = $fname"
    echo "compiler        = $compiler"
    echo "bin_dir         = $bin_dir"
    echo "nrepeat         = $nrepeat"
    echo "niter           = $niter"
    echo "nx              = $nx"
    echo "ny              = $ny"
    echo "nz              = $nz"
    echo "lx              = $lx"
    echo "ly              = $ly"
    echo "lz              = $lz"
    echo "min_ntasks      = $min_ntasks"
    echo "inc_ntasks      = $inc_ntasks"
    echo "max_ntasks      = $max_ntasks"
    if ! test "$subcomm" = "true"; then
        subcomm="false"
    fi
    echo "subcomm         = $subcomm"
    echo "enable_caf      = $enable_caf"
    echo "--------------------------------"

    mkdir -p -v "$compiler"
    cd "$compiler"

    ntasks=$min_ntasks
    while (($ntasks <= $max_ntasks)); do
        nodes=$((ntasks/ntasks_per_node))

        # avoid nodes = 0
        if test $nodes = 0; then
            nodes=1
        fi

        echo "Submit job with $ntasks tasks on $nodes nodes using the $compiler version"

        if test "$enable_caf" = "yes"; then
            fn="submit_caf_random_nx_${nx}_ny_${ny}_nz_${nz}_nodes_${nodes}.sh"
        else
            fn="submit_random_nx_${nx}_ny_${ny}_nz_${nz}_nodes_${nodes}.sh"
        fi

        cp "../$fname" $fn
        sed -i "s:JOBNAME:$compiler-random:g" $fn
        sed -i "s:COMPILER:$compiler:g" $fn

        sed -i "s:NREPEAT:$nrepeat:g" $fn
        sed -i "s:NODES:$nodes:g" $fn
        sed -i "s:NTASKS:$ntasks:g" $fn
        sed -i "s:--niter NITER:--niter $niter:g" $fn
        sed -i "s:NX:$nx:g" $fn
        sed -i "s:NY:$ny:g" $fn
        sed -i "s:NZ:$nz:g" $fn

        sed -i "s:--lx LX:--lx $lx:g" $fn
        sed -i "s:--ly LY:--ly $ly:g" $fn
        sed -i "s:--lz LZ:--lz $lz:g" $fn
        sed -i "s:--xlen LX:--xlen $lx:g" $fn
        sed -i "s:--ylen LY:--ylen $ly:g" $fn
        sed -i "s:--zlen LZ:--zlen $lz:g" $fn

        sed -i "s:BIN_DIR:$bin_dir:g" $fn
        sed -i "s:ENABLE_CAF:$enable_caf:g" $fn
        sed -i "s:SUBCOMM:$subcomm:g" $fn

        sbatch $fn

        ntasks=$((ntasks*2))
    done

    cd ..
}

# Argument order of run_jobs
# machine
# ntasks_per_node
# compiler : cray or gnu
# bin_dir
# nrepeat
# niter
# nx
# ny
# nz
# lx
# ly
# lz
# min_ntasks
# inc_ntasks
# max_ntasks
# subcomm
# enable_caf : "no" or "yes"

print_help() {
    echo "Script to submit strong / weak scaling jobs"
    echo "where the number of cores is doubled in each"
    echo "iteration from '-l' to '-u'"
    echo "Arguments:"
    echo "    -m    machine to run on, either 'cirrus' or 'archer2'"
    echo "    -h    print this help message"
    echo "    -l    lower bound of cores"
    echo "    -j    increment of cores"
    echo "    -u    upper bound of cores"
    echo "    -r    number of repetitions"
    echo "    -i    number of iterations per repetition"
    echo "    -x    number of grid cells in the horizontal direction x"
    echo "    -y    number of grid cells in the horizontal direction y"
    echo "    -z    number of grid cells in the vertical direction z"
    echo "    -a    domain extent in the horizontal direction x"
    echo "    -b    domain extent in the horizontal direction y"
    echo "    -c    domain extent in the vertical direction z"
    echo "    -s    use sub-communicator (optional)"
}

machine=''
subcomm="false"
while getopts "h?m:l:u:j:r:i:x:y:z:a:b:c:s": option; do
    case "$option" in
        a)
            lx=$OPTARG
            ;;
        b)
            ly=$OPTARG
            ;;
        c)
            lz=$OPTARG
            ;;
        h|\?)
            print_help
            exit 0
            ;;
        i)
            niter=$OPTARG
            ;;
        j)
            inc_cores=$OPTARG
            ;;
        l)
            min_cores=$OPTARG
            ;;
    	m)
            machine=$OPTARG
      	    ;;
        r)
            nrep=$OPTARG
            ;;
        s)
            subcomm="true"
            ;;
        u)
            max_cores=$OPTARG
            ;;
        x)
            nx=$OPTARG
            ;;
        y)
            ny=$OPTARG
            ;;
        z)
            nz=$OPTARG
            ;;
    esac
done

if ! test "$machine" = "archer2" && ! test "$machine" = "cirrus"; then
    echo "Only 'archer2' and 'cirrus' machines supported. Exiting."
    exit 1
fi

# set bin directories
source "../$machine.sh"

echo "Submiting jobs on $machine with $min_cores to $max_cores cores."
echo "Each job is repeated $nrep times with $niter iterations per repetition."

j=0
for bin_dir in ${bins[*]}; do
    compiler="${compilers[$j]}"
    with_caf="${enable_caf[$j]}"

    run_jobs $machine $ntasks_per_node $compiler "$bin_dir" $nrep $niter $nx $ny $nz $lx $ly $lz $min_cores $inc_cores $max_cores $subcomm $with_caf

    j=$((j+1))
done
