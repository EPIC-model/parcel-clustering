#!/bin/bash

run_jobs() {

    local machine=${1}
    local fname="submit_${machine}_read.sh"
    local ntasks_per_node=${2}

    local compiler=${3}

    local bin_dir=${4}
    local nrepeat=${5}
    local niter=${6}
    local fullname=${7}
    local offset=${8}
    local nfiles=${9}

    local min_ntasks=${10}
    local inc_ntasks=${11}
    local max_ntasks=${12}
    local subcomm=${13}
    local size_factor=${14}
    local time_limit=${15}

    local bname=$(basename $fullname)
    local dname=$(dirname $fullname)
    if [[ $bname =~ ^[a-z_]*([0-9]*)x([0-9]*)x([0-9]*)[_a-z]* ]]; then
        local nx=${BASH_REMATCH[1]}
        local ny=${BASH_REMATCH[2]}
        local nz=${BASH_REMATCH[3]}
        echo "$nx, $ny, $nz"
    else
        echo "Error in matching: $bname"
        exit 1
    fi

    if [[ $bname == *"early"* ]]; then
        name_tag="early"
    elif [[ $bname == *"late"* ]]; then
        name_tag="late"
    else
        echo "Error: Could not identify if 'early' or 'late' simulation time. Exiting."
        exit 1
    fi

    echo "--------------------------------"
    echo "Run jobs with following options:"
    echo "machine         = $machine"
    echo "ntasks_per_node = $ntasks_per_node"
    echo "fname           = $fname"
    echo "compiler        = $compiler"
    echo "bin_dir         = $bin_dir"
    echo "nrepeat         = $nrepeat"
    echo "niter           = $niter"
    echo "dirname         = $dname"
    echo "basename        = $bname"
    echo "name_tag        = $name_tag"
    echo "size_factor     = $size_factor"
    echo "nx              = $nx"
    echo "ny              = $ny"
    echo "nz              = $nz"
    echo "offset          = $offset"
    echo "nfiles          = $nfiles"
    echo "min_ntasks      = $min_ntasks"
    echo "inc_ntasks      = $inc_ntasks"
    echo "max_ntasks      = $max_ntasks"
    echo "time_limit      = $time_limit"
    if ! test "$subcomm" = "true"; then
        subcomm="false"
    fi
    echo "subcomm         = $subcomm"
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

        fn="submit_${machine}_read_${name_tag}_nx_${nx}_ny_${ny}_nz_${nz}_nodes_${nodes}.sh"

        cp "../$fname" $fn
        sed -i "s:JOBNAME:$machine-$compiler-read:g" $fn
        sed -i "s:COMPILER:$compiler:g" $fn
        sed -i "s:MACHINE:$machine:g" $fn
        sed -i "s:NAMETAG:$name_tag:g" $fn
        sed -i "s/--time=TIMELIMIT/--time=$time_limit/g" $fn

        sed -i "s:NX:$nx:g" $fn
        sed -i "s:NY:$ny:g" $fn
        sed -i "s:NZ:$nz:g" $fn

        sed -i "s:NREPEAT:$nrepeat:g" $fn
        sed -i "s:NODES:$nodes:g" $fn
        sed -i "s:--ntasks=NTASKS:--ntasks=$ntasks:g" $fn
        sed -i "s:-np NTASKS:-np $ntasks:g" $fn
        sed -i "s:--niter NITER:--niter $niter:g" $fn
        sed -i "s:--dirname DIRNAME:--dirname $dname:g" $fn
        sed -i "s:--ncbasename NC_BASENAME:--ncbasename $bname:g" $fn
        sed -i "s:--offset OFFSET:--offset $offset:g" $fn
        sed -i "s:--nfiles NFILES:--nfiles $nfiles:g" $fn
        sed -i "s:--size-factor SIZE_FACTOR:--size-factor $size_factor:g" $fn

        sed -i "s:BIN_DIR:$bin_dir:g" $fn
        sed -i "s:SUBCOMM:$subcomm:g" $fn

        sbatch $fn

        ntasks=$((ntasks*inc_ntasks))
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
# basename
# offset
# nfiles
# min_ntasks
# inc_ntasks
# max_ntasks
# subcomm
# size_factor

print_help() {
    echo "Script to submit strong scaling jobs reading netCDF files"
    echo "where the number of cores is doubled in each"
    echo "iteration from '-l' to '-u'"
    echo ""
    echo "Requirement:"
    echo "          EPIC parcel files (netCDF format)."
    echo ""
    echo "Arguments:"
    echo "    -m    machine to run on, e.g. 'cirrus', 'archer2', 'hotlum'"
    echo "          (requirement: <machine>.sh and 'submit_<machine>_read.sh)"
    echo "    -h    print this help message"
    echo "    -l    lower bound of cores"
    echo "    -j    increment of cores"
    echo "    -u    upper bound of cores"
    echo "    -r    number of repetitions"
    echo "    -i    number of iterations per repetition"
    echo "    -b    basename of EPIC parcel netCDF files,"
    echo "          e.g. <epic_sim> for <epic_sim>_0000000011_parcels.nc"
    echo "    -o    starting index of EPIC parcel files, e.g."
    echo "          -o 11 for <epic_sim>_0000000011_parcels.nc"
    echo "    -n    number of EPIC netCDF files;"
    echo "          if the number of iterations > number of EPIC files, the"
    echo "          file reading cycle restarts at the starting index"
    echo "    -s    use sub-communicator (optional)"
    echo "    -f    parcel container size factor, --size-factor"
    echo "    -t    time limit of job, default: 00:30:00"
}

check_for_input() {
    if ! test "${2}"; then
        echo "Please specify '${1}'. Exiting."
	exit 1
    fi
}


machine=''

# default options:
subcomm="false"
inc_cores=2
nrep=1
niter=1
timelimit="00:30:00"

while getopts "h?m:l:u:j:r:i:b:o:n:sf:t:" option; do
    case "$option" in
        b)
            file_base_name=$OPTARG
            ;;
        h|\?)
            print_help
            exit 0
            ;;
        f)
            szf=$OPTARG
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
        n)
            num_files=$OPTARG
            ;;
    	m)
            machine=$OPTARG
      	    ;;
        o)
            file_offset=$OPTARG
            ;;
        r)
            nrep=$OPTARG
            ;;
        s)
            subcomm="true"
            ;;
        t)
            timelimit=$OPTARG
            ;;
        u)
            max_cores=$OPTARG
            ;;
    esac
done

if ! test -f "../$machine.sh"; then
    echo "Unable to run on $machine. The file ${machine}.sh does not exist. Exiting."
    exit 1
fi

if ! test -f "submit_${machine}_read.sh" ; then
    echo "Unable to run on $machine. The file submit_${machine}_read.sh does not exist. Exiting."
    exit 1
fi

# set bin directories
source "../$machine.sh"

check_for_input "ntasks_per_node" $ntasks_per_node
check_for_input "nrep" $nrep
check_for_input "niter" $niter
check_for_input "file_base_name" $file_base_name
check_for_input "file_offset" $file_offset
check_for_input "num_files" $num_files
check_for_input "min_cores" $min_cores
check_for_input "inc_cores" $inc_cores
check_for_input "max_cores" $max_cores
check_for_input "subcomm" $subcomm
check_for_input "size_factor" $szf
check_for_input "time_limit" $timelimit

echo "Submiting jobs on $machine with $min_cores to $max_cores cores."
echo "Each job is repeated $nrep times with $niter iterations per repetition."

j=0
for bin_dir in ${bins[*]}; do
    compiler="${compilers[$j]}"

    run_jobs $machine $ntasks_per_node $compiler "$bin_dir" $nrep $niter $file_base_name $file_offset $num_files $min_cores $inc_cores $max_cores $subcomm $szf $timelimit

    j=$((j+1))
done
