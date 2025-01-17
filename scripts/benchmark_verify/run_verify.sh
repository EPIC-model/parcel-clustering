#!/bin/bash

fname="submit_verify.sh"

for c in "gnu" "cray"; do
    mkdir "$c"
    cd "$c"
    for i in "p2p" "rma" "shmem"; do
	mkdir "$i"
	cd "$i"
	cp "../../$fname" .

	sed -i "s:#SBATCH --job-name=JOBNAME:#SBATCH --job-name=$c-$i:g" $fname
	sed -i "s:COMPILER:$c:g" $fname
	sed -i "s:GRAPH_TYPE:$i:g" $fname
	
	echo "Submit job $i with $c"
	sbatch $fname
	cd ..
    done
    cd ..
done
