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
instructios of the [Cirrus MPI build](https://github.com/hpc-uk/build-instructions/blob/main/libs/openmpi/build_openmpi_5.0.0_cirrus_gcc10.md).

```bash
module load libtool/2.4.7
module load gcc/10.2.0
export MPI_ROOT=/work/e710/e710/mf248/gcc/10.2.0/openmpi/5.0.7
export PATH=$MPI_ROOT/bin:$PATH
export MPIF90=$MPI_ROOT/bin/mpif90
export MPIFORT=$MPI_ROOT/bin/mpifort
export MPICC=$MPI_ROOT/bin/mpicc
export MPICXX=$MPI_ROOT/bin/mpic++
export NETCDF_C_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
export NETCDF_FORTRAN_DIR=/work/e710/e710/mf248/gcc/10.2.0/netcdf
export PATH=$PATH:$NETCDF_C_DIR/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$NETCDF_C_DIR/lib:$MPI_ROOT/lib
export CPLUS_INCLUDE_PATH=$NETCDF_C_DIR/include:$MPI_ROOT/include:$CPLUS_INCLUDE_PATH
export C_INCLUDE_PATH=$NETCDF_C_DIR/include:$MPI_ROOT/include:$C_INCLUDE_PATH
CC=$MPICC CXX=$MPICXX FC=$MPIFORT ../configure --prefix=/work/e710/e710/mf248/gnu
```

UCX version: 1.16.0
HWLOC version: 2.9.3
PMIX version: 5.0.6
OpenMPI version: 5.0.7
```bash
Open MPI configuration:
-----------------------
Version: 5.0.7
MPI Standard Version: 3.1
Build MPI C bindings: yes
Build MPI Fortran bindings: mpif.h, use mpi, use mpi_f08
Build MPI Java bindings (experimental): no
Build Open SHMEM support: yes
Debug build: no
Platform file: (none)
 
Miscellaneous
-----------------------
Atomics: GCC built-in style atomics
Fault Tolerance support: mpi
HTML docs and man pages: installing packaged docs
hwloc: external
libevent: external
Open UCC: no
pmix: external
PRRTE: internal
Threading Package: pthreads
 
Transports
-----------------------
Cisco usNIC: no
Cray uGNI (Gemini/Aries): no
Intel Omnipath (PSM2): no (not found)
Open UCX: yes
OpenFabrics OFI Libfabric: no (not found)
Portals4: no (not found)
Shared memory/copy in+copy out: yes
Shared memory/Linux CMA: yes
Shared memory/Linux KNEM: yes
Shared memory/XPMEM: no
TCP: yes
 
Accelerators
-----------------------
CUDA support: no
ROCm support: no
 
OMPIO File Systems
-----------------------
DDN Infinite Memory Engine: no
Generic Unix FS: yes
IBM Spectrum Scale/GPFS: no (not found)
Lustre: no (not found)
PVFS2/OrangeFS: no
```
