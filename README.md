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

### Installation on ARCHER2
Please also consult the [ARCHER2 documentation](https://docs.archer2.ac.uk).

#### Cray Compiling Environment (CCE) suite
```bash
module load cce/15.0.0
module load cray-mpich/8.1.23
module load cray-hdf5-parallel/1.12.2.1
module load cray-openshmemx/11.5.7
module load cray-netcdf-hdf5parallel/4.9.0.1
module load cpe/23.09
export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_DIR/crayclang/14.0
export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_DIR/crayclang/14.0
export CXX=CC
export CC=cc
export FC=ftn
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
```

#### GNU Compiler Collection (GCC) suite
```bash
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load cray-openshmemx
module load PrgEnv-gnu
module load cray-hdf5-parallel
module load cray-netcdf-hdf5parallel
module load cray-openshmemx
module load load-epcc-module;
module load  extra-compilers/1.0
module load cpe/23.09
export NETCDF_C_DIR=$NETCDF_DIR
export NETCDF_FORTRAN_DIR=$NETCDF_DIR
export CXX=CC
export CC=cc
export FC=ftn
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
```


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
CXX="mpicxx -cxx=icpc -lsma" CC="mpicc -cc=icc -lsma" FC="mpif08 -fc=ifort -lsma" ../configure --prefix=$PREFIX
```

<!-- #### Intel MPI with Intel compiler suite
```bash
module load intel-20.4/mpi
module load intel-20.4/compilers
module load netcdf-parallel/4.9.2-intel20-impi20
export NETCDF_C_DIR=$NETCDF_DIR
export NETCDF_FORTRAN_DIR=$NETCDF_DIR
CXX=mpiicpc CC=mpiicc FC=mpiifort ../configure
```
-->

#### OpenMPI with GNU compiler suite

That is the build with OpenMPI/4.1.6.
```bash
module load libtool/2.4.7
module load gcc/10.2.0
module load openmpi/4.1.6
module load hdf5parallel/1.14.3-gcc10-ompi416
export NETCDF_C_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
export NETCDF_FORTRAN_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
export PATH=$PATH:$NETCDF_C_DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_C_DIR/lib
export CPLUS_INCLUDE_PATH=$NETCDF_C_DIR/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=$NETCDF_C_DIR/include:$C_INCLUDE_PATH
CC=mpicc CXX=mpic++ FC=mpifort ../configure --prefix=/work/e710/e710/mf248/gnu
```

However, we use the latest version OpenMPI/5.0.7 which we build following the
instructios of the [Cirrus documentation](https://github.com/hpc-uk/build-instructions/blob/main/libs/openmpi/build_openmpi_4.1.6_cirrus_gcc10.md).
