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


