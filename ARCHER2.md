## Installation on ARCHER2
Please also consult the [ARCHER2 documentation](https://docs.archer2.ac.uk).

### Cray Compiling Environment (CCE) suite
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

### GNU Compiler Collection (GCC) suite
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

