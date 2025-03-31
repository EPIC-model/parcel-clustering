#!/bin/bash

#gnu_bin=/lus/bnchlu1/shanks/EPIC/parcel-clustering/build_gnu/bin
cray_bin=/lus/bnchlu1/shanks/EPIC/parcel-clustering/build_cce/bin

cray_compiler="cray"
#gnu_compiler="gnu"

declare -a bins=("$cray_bin") #"$gnu_bin"
declare -a compilers=("$cray_compiler") #"$gnu_compiler"

ntasks_per_node=128
