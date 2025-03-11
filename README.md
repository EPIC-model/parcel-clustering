# Parallel parcel / particle cluster algorithm

## Requirements
* GNU Autotools
    - GNU Autoconf >= 2.69
    - GNU automake >= 1.16.1
    - GNU libtool >= 2.4.6
* GNU Fortran >= 9.3.0
* netCDF >= 4.8.1
* netCDF-Fortran >= 4.5.4
* HDF5 >= 1.12.1
* MPI (MPICH, Cray MPI, OpenMPI)
* OpenSHMEM (included in OpenMPI)

## Installation
1. `./bootstrap`
2. `mkdir build; cd build`
3. `CXX=mpic++ CC=mpicc FC=mpifort ../configure --prefix=$PREFIX`
4. `make; make install`

Note: The `$PREFIX` variable denotes the installation directory. The flags `CXX`, `CC` and `FC` may also be other compiler wrappers.

* [Installation on ARCHER2](ARCHER2.md)
* [Installation on Cirrus](Cirrus.md)

## How to run the verification benchmark
In order to run the code verification benchmark you need to add two machine-specific files
to the directory [scripts/benchmark_verify](https://github.com/EPIC-model/parcel-clustering/tree/update-doc/scripts/benchmark_verify) directory:

* `<machine>.sh`
* `submit_<machine>_verify.sh`

where `<machine>` is the name of the computing system (or any arbitrary name). As an example of a submission script see [submit_archer2_verify.sh](scripts/benchmark_verify/submit_archer2_verify.sh).
For the verification benchmark, the submission script must specify the following placeholders:

| placeholder | description                               |
| ----------- | ----------------------------------------- |
| COMPILER    | name of the compiler suite, e.g. cce, gnu |
| COMM_TYPE   | communication layer, e.g. shmem, p2p, rma |
| N_SAMPLES   | number of random samples                  |
| SEED        | seed for the random sample generator      |
| BIN_DIR     | bin directory of the executable           |

The file `<machine.sh>` must contain the variable `ntasks_per_node` and two arrays
`bins` and `compilers` that specify the location of the executables and the name of the compiler suite, respectively.
As an example consult [archer2.sh](scripts/archer2.sh). Once these files are specified the verification benchmark can be
started within the directory [scripts/benchmark_verify](scripts/benchmark_verify) using the following command
```bash
bash run_verify.sh -s [seed] -n [number of samples] -m [machine] -c [communication layer]
```
where the tags within square brackets must be specified.
