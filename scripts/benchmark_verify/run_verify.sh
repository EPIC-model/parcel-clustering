#!/bin/bash

fname="submit_verify.sh"

n_samples=1000
seed=42

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
        sed -i "s:N_SAMPLES:$n_samples:g" $fname
        sed -i "s:SEED:$seed:g" $fname

        echo "Submit job $i with $c. Running $n_samples samples with seed $seed."
        sbatch $fname
        cd ..
    done
    cd ..
done

# Run Coarray Fortran
cd "cray"
mkdir "caf"
cd "caf"
cp "../../$fname" .

sed -i "s:#SBATCH --job-name=JOBNAME:#SBATCH --job-name=cray-caf:g" $fname
sed -i "s:COMPILER:cray-caf:g" $fname
sed -i "s:GRAPH_TYPE:caf:g" $fname
sed -i "s:N_SAMPLES:$n_samples:g" $fname
sed -i "s:SEED:$seed:g" $fname

echo "Submit job caf with cray. Running $n_samples samples with seed $seed."
sbatch $fname
cd ..
