SUBDIRS = mpi netcdf utils parcels

AM_FCFLAGS =				\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf	\
	-I $(top_builddir)/src/parcels

bin_PROGRAMS = benchmark_random	\
               benchmark_read

if ENABLE_PYTHON
bin_PROGRAMS += benchmark_verify
endif

benchmark_random_SOURCES = utils.f90 benchmark_random.f90

benchmark_random_LDADD =				\
	$(top_builddir)/src/utils/libepic_utils.la      \
        $(top_builddir)/src/mpi/libepic_mpi.la          \
        $(top_builddir)/src/netcdf/libepic_netcdf.la	\
        $(top_builddir)/src/parcels/libparcels.la

benchmark_read_SOURCES = utils.f90 benchmark_read.f90

benchmark_read_LDADD = 					\
	$(top_builddir)/src/utils/libepic_utils.la      \
        $(top_builddir)/src/mpi/libepic_mpi.la          \
        $(top_builddir)/src/netcdf/libepic_netcdf.la	\
        $(top_builddir)/src/parcels/libparcels.la


if ENABLE_PYTHON
benchmark_verify_SOURCES =			\
	utils.f90	           		\
	serial/parcel_nearest_serial.f90	\
	serial/parcel_merge_serial.f90         	\
	benchmark_verify.f90

benchmark_verify_LDADD =				\
	$(top_builddir)/src/utils/libepic_utils.la      \
        $(top_builddir)/src/mpi/libepic_mpi.la          \
        $(top_builddir)/src/netcdf/libepic_netcdf.la	\
        $(top_builddir)/src/parcels/libparcels.la
endif

clean-local:
	rm -f *.mod *.smod
