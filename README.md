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
3. `../configure --prefix=$PREFIX --enable-python`
4. `make; make install`

Note: The `$PREFIX` variable denotes the installation directory.

### Installation on Cirrus
Please also consult the [Cirrus documentation](https://docs.cirrus.ac.uk/user-guide/development/).

#### HPE MPT MPI with Intel compiler suite
```bash
module load libtool/2.4.7
module load mpt/2.25
module load intel-20.4/compilers
module load netcdf-parallel/4.9.2-intel20-mpt225
export NETCDF_FORTRAN_DIR=$NETCDF_DIR
export NETCDF_C_DIR=$NETCDF_DIR
export MPICC_CC=icc
export MPICXX_CXX=icpc
CXX="mpicxx -cxx=icpc -lsma" CC="mpicc -cc=icc -lsma" FC="mpif08 -fc=ifort -lsma" ../configure --prefix=$PREFIX --enable-python
```

#### Intel MPI with Intel compiler suite
```bash
module load intel-20.4/mpi
module load intel-20.4/compilers
module load netcdf-parallel/4.9.2-intel20-impi20
export NETCDF_C_DIR=$NETCDF_DIR
export NETCDF_FORTRAN_DIR=$NETCDF_DIR
CXX=mpiicpc CC=mpiicc FC=mpiifort ../configure --enable-python
```

#### OpenMPI with GNU compiler suite
```bash
module load libtool/2.4.7
module load gcc/10.2.0
module load openmpi/4.1.6
module load hdf5parallel/1.14.3-gcc10-ompi416
export NETCDF_C_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
export NETCDF_FORTRAN_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
export PATH=$PATH:$NETCDF_C_DIR/bin
export MPI_DIR=/work/y07/shared/cirrus-software/openmpi/4.1.6
CC=mpicc CXX=mpicxx FC=mpif90 ../configure --enable-python
```
