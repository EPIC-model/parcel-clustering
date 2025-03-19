#!/bin/bash

print_help() {
    echo "Script to prepare RT test case."
    echo "Arguments:"
    echo "    -h    print this help message"
    echo "    -n    number of comput nodes to run on"
    echo "    -x    number of grid cells in the horizontal direction x"
    echo "    -y    number of grid cells in the horizontal direction y"
    echo "    -z    number of grid cells in the vertical direction z"
    echo "    -e    early time limit; for early restart run"
    echo "    -l    late time limit; for late restart run"
}

check_for_input() {
    if ! test "${2}"; then
        echo "Please specify '${1}'. Exiting."
        exit 1
    fi
}

while getopts "h?x:y:z:n:e:l:": option; do
    case "$option" in
      h|\?)
            print_help
            exit 0
            ;;
        e)
            early_limit=$OPTARG
            ;;
        l)
            late_limit=$OPTARG
            ;;
        n)
            nodes=$OPTARG
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

check_for_input "nx" $nx
check_for_input "ny" $ny
check_for_input "nz" $nz
check_for_input "nodes" $nodes
check_for_input "early limit" $early_limit
check_for_input "late limit" $late_limit

ntasks=$((128*$nodes))

echo "Prepare RT simulation with mesh ${nx}x${ny}x${nz}."

python rayleigh_taylor.py --nx ${nx} --ny ${ny} --nz ${nz} --epsilon 0.1 --ape-calculation "none"

sim_dir="rt-${nx}x${ny}x${nz}"
mkdir -p "${sim_dir}"
cp rt-template/* "${sim_dir}/"

cd ${sim_dir}
mv "../rt_${nx}x${ny}x${nz}.nc" .
mkdir "early-time"
mkdir "late-time"
sed -i "s:NXxNYxNZ:${nx}x${ny}x${nz}:g" rt_epic_prepare.config
sed -i "s:NXxNYxNZ:${nx}x${ny}x${nz}:g" rt_epic_early_restart.config
sed -i "s:NXxNYxNZ:${nx}x${ny}x${nz}:g" rt_epic_late_restart.config

sed -i "s:EARLY_LIMIT:${early_limit}:g" rt_epic_early_restart.config
sed -i "s:LATE_LIMIT:${late_limit}:g" rt_epic_late_restart.config

sed -i "s:NODES:${nodes}:g" submit_rt_epic_prepare.sh
sed -i "s:NTASKS:${ntasks}:g" submit_rt_epic_prepare.sh

sed -i "s:NODES:${nodes}:g" submit_rt_epic_early_time_restart.sh
sed -i "s:NTASKS:${ntasks}:g" submit_rt_epic_early_time_restart.sh
sed -i "s:NXxNYxNZ:${nx}x${ny}x${nz}:g" submit_rt_epic_early_time_restart.sh

sed -i "s:NODES:${nodes}:g" submit_rt_epic_late_time_restart.sh
sed -i "s:NTASKS:${ntasks}:g" submit_rt_epic_late_time_restart.sh
sed -i "s:NXxNYxNZ:${nx}x${ny}x${nz}:g" submit_rt_epic_late_time_restart.sh
cd ..
