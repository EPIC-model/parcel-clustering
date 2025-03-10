#!/bin/bash

gnu_bin="/work/e710/e710/mf248/gnu/bin"
cray_bin="/work/e710/e710/mf248/cray/bin"

cray_compiler="cray"
gnu_compiler="gnu"

declare -a bins=("$gnu_bin" "$cray_bin")
declare -a compilers=("$gnu_compiler" "$cray_compiler")

ntasks_per_node=128
