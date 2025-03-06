#!/bin/bash

gnu_bin="/work/e710/e710/mf248/gnu/bin"
intel_mpi_bin="/work/e710/e710/mf248/intel-mpi/bin"
intel_mpt_bin="/work/e710/e710/mf248/intel-mpt/bin"

declare -a bins=("$gnu_bin" "$intel_mpi_bin" "$intel_mpt_bin")
declare -a compilers=("gnu" "intel-mpi" "intel-mpt")
declare -a enable_caf=("no" "no" "no")

ntasks_per_node=36
