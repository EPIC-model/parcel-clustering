## Installation on Hotlum

### Cray Compiling Environment (CCE) suite
```bash
module load PrgEnv-cray
module load cray-hdf5-parallel
module load cray-openshmemx
module load cray-netcdf-hdf5parallel
module load -f cpe/24.11
export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
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
module load -f cpe/24.11
export NETCDF_C_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export NETCDF_FORTRAN_DIR=$CRAY_NETCDF_HDF5PARALLEL_PREFIX
export CXX=CC
export CC=cc
export FC=ftn
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH
```

