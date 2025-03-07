#!/bin/bash

gnu_bin="/work/e710/e710/mf248/gnu/clustering/bin"
cray_bin="/work/e710/e710/mf248/cray/clustering/bin"
caf_bin="/work/e710/e710/mf248/cray-caf/clustering/bin"

cray_compiler="cray"
gnu_compiler="gnu"
caf_compiler=$cray_compiler

declare -a bins=("$gnu_bin" "$cray_bin" "$caf_bin")
declare -a compilers=("$gnu_compiler" "$cray_compiler" "$caf_compiler")
declare -a enable_caf=("no" "no" "yes")

ntasks_per_node=128
