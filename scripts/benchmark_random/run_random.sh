#!/bin/bash

run_jobs() {

    # 22 Jan 2025
    # https://stackoverflow.com/a/6212408
    # https://stackoverflow.com/a/12128447


    local fname=${1}

    # "gnu" or "cray"
    local compiler=${2}

    mkdir -p $compiler; cd $compiler

    local subcomm="false"

    local nrepeat=${3}
    local niter=${4}
    local nx=${5}
    local ny=${6}
    local lx=${7}
    local ly=${8}
    local begin=${9}
    local end=${10}
    local subcomm=${11}

    echo "--------------------------------"
    echo "Run jobs with following options:"
    echo "fname    = $fname"
    echo "compiler = $compiler"
    echo "nrepeat  = $nrepeat"
    echo "niter    = $niter"
    echo "nx       = $nx"
    echo "ny       = $ny"
    echo "lx       = $lx"
    echo "ly       = $ly"
    echo "begin    = $begin"
    echo "end      = $end"
    if ! test "$subcomm" = "true"; then
        subcomm="false"
    fi
    echo "subcomm  = $subcomm"
    echo "--------------------------------"

    for i in $(seq $begin 1 $end); do
        ntasks=$((2**i))
        nodes=$((ntasks/128))
    
        # 2 March 2024
        # https://stackoverflow.com/a/10415158
        nodes=$((nodes==0 ? 1 : nodes))
        echo "Submit job with $ntasks tasks on $nodes nodes using the $compiler version"
    
        fn="submit_random_nx_${nx}_ny_${ny}_nodes_${nodes}.sh"
        cp ../$fname $fn
        sed -i "s:JOBNAME:$compiler-random:g" $fn
        sed -i "s:COMPILER:$compiler:g" $fn

        sed -i "s:NREPEAT:$nrepeat:g" $fn
        sed -i "s:NODES:$nodes:g" $fn
        sed -i "s:NTASKS:$ntasks:g" $fn
        sed -i "s:--niter NITER:--niter $niter:g" $fn
        sed -i "s:NX:$nx:g" $fn
        sed -i "s:NY:$ny:g" $fn
        sed -i "s:--lx LX:--lx $lx:g" $fn
        sed -i "s:--ly LY:--ly $ly:g" $fn
        sed -i "s:--xlen LX:--xlen $lx:g" $fn
        sed -i "s:--ylen LY:--ylen $ly:g" $fn

        sed -i "s:SUBCOMM:$subcomm:g" $fn
    
        sbatch $fn
    done
    
    cd ..
}

# Argument order:
# fname
# compiler : cray or gnu
# nrepeat
# niter
# nx
# ny
# lx
# ly
# begin
# end
# subcomm

for compiler in "cray"; do #"gnu" "cray"; do
    # 1 node to 8 nodes
    run_jobs "submit_random.sh" $compiler 1 5 256 512 80 160 7 10 "false"

    # 2 nodes to 32 nodes
    run_jobs "submit_random.sh" $compiler 1 5 512 512 160 160 8 12 "false"

    # 8 nodes to 128 nodes
    run_jobs "submit_random.sh" $compiler 1 5 1024 1024 320 320 10 14 "false"
done
