#!/bin/bash

run_jobs() {

    local machine=${1}
    local fname="submit_${machine}_random.sh"

    # "gnu" or "cray"
    local compiler=${2}

    local bin_dir=${3}
    local nrepeat=${4}
    local niter=${5}
    local nx=${6}
    local ny=${7}
    local nz=${8}
    local lx=${9}
    local ly=${10}
    local lz=${11}
    local begin=${12}
    local end=${13}
    local subcomm=${14}
    local enable_caf=${15}

    echo "--------------------------------"
    echo "Run jobs with following options:"
    echo "machine    = $machine"
    echo "fname      = $fname"
    echo "compiler   = $compiler"
    echo "bin_dir    = $bin_dir"
    echo "nrepeat    = $nrepeat"
    echo "niter      = $niter"
    echo "nx         = $nx"
    echo "ny         = $ny"
    echo "nz         = $nz"
    echo "lx         = $lx"
    echo "ly         = $ly"
    echo "lz         = $lz"
    echo "begin      = $begin"
    echo "end        = $end"
    if ! test "$subcomm" = "true"; then
        subcomm="false"
    fi
    echo "subcomm    = $subcomm"
    echo "enable_caf = $enable_caf"
    echo "--------------------------------"

    mkdir -p -v "$compiler"
    cd "$compiler"

    for i in $(seq $begin 1 $end); do
        ntasks=$((2**i))
        nodes=$((ntasks/128))

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
    done

    cd ..
}

# Argument order:
# fname
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
# begin
# end
# subcomm
# enable_caf : "no" or "yes"

machine="archer2"
gnu_bin="/work/e710/e710/mf248/gnu/clustering/bin"
cray_bin="/work/e710/e710/mf248/cray/clustering/bin"
caf_bin="/work/e710/e710/mf248/cray-caf/clustering/bin"

cray_compiler="cray"
gnu_compiler="gnu"
caf_compiler=$cray_compiler

compiler=$gnu_compiler
enable_caf="no"

for bin_dir in $gnu_bin $cray_bin $caf_bin; do
    if ! test -d "$bin_dir"; then
        echo "No bin directory: $bin_dir"
        break
    fi

    if test "$bin_dir" = "$gnu_bin"; then
        compiler=$gnu_compiler
        enable_caf="no"
    fi

    if test "$bin_dir" = "$cray_bin"; then
        compiler=$cray_compiler
        enable_caf="no"
    fi

    if test "$bin_dir" = "$caf_bin"; then
        compiler=$caf_compiler
        enable_caf="yes"
    fi

    run_jobs "$machine" "$compiler" "$bin_dir" 1 5 256 512 64 80 160 20 7 10 "false" "$enable_caf"

    # 2 nodes to 32 nodes
    run_jobs "$machine" "$compiler" "$bin_dir" 1 5 512 512 64 160 160 20 8 12 "false" "$enable_caf"

    # 8 nodes to 128 nodes
    run_jobs "$machine" "$compiler" "$bin_dir" 1 5 1024 1024 64 320 320 20 10 14 "false" "$enable_caf"
    
done
