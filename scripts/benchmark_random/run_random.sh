#!/bin/bash

fname="submit_random.sh"

# "gnu" or "cray"
compiler="cray"

mkdir -p $compiler; cd $compiler

nx=1024
ny=1024
lx=320
ly=320

for i in $(seq 7 1 14); do
	ntasks=$((2**i))
	nodes=$((ntasks/128))

	# 2 March 2024
    # https://stackoverflow.com/a/10415158
    nodes=$((nodes==0 ? 1 : nodes))
    echo "Submit job with $ntasks tasks on $nodes nodes using the $compiler version"

	fn="submit_random_$nodes.sh"
	cp ../$fname $fn
	sed -i "s:JOBNAME:$compiler-random:g" $fn
	sed -i "s:COMPILER:$compiler:g" $fn
	
	sed -i "s:NODES:$nodes:g" $fn
	sed -i "s:NTASKS:$ntasks:g" $fn
	sed -i "s:--nx NX:--nx $nx:g" $fn
    sed -i "s:--ny NY:--ny $ny:g" $fn
    sed -i "s:--lx LX:--lx $lx:g" $fn
    sed -i "s:--ly LY:--ly $ly:g" $fn
    sed -i "s:--xlen LX:--xlen $lx:g" $fn
    sed -i "s:--ylen LY:--ylen $ly:g" $fn

	sbatch $fn
done

cd ..
