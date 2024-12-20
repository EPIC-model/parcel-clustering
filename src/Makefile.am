SUBDIRS = mpi netcdf utils

AM_FCFLAGS =				\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf

bin_PROGRAMS = benchmark_parcel_merging
benchmark_parcel_merging_SOURCES =		\
	parcels/parcel_ops.f90			\
	parcels/parcel_container.f90 		\
	parcels/parcel_ellipsoid.f90		\
	parcels/parcels.f90			\
	parcels/parcel_bc.f90                   \
	parcels/parcel_mpi.f90			\
	parcels/parcel_split.f90		\
	parcels/parcel_init.f90			\
	parcels/parcel_nearest_graph.f90	\
	parcels/parcel_nearest.f90		\
	parcels/parcel_merge.f90		\
	parcels/parcel_netcdf.f90		\
	utils.f90				\
	benchmark_parcel_merging.f90

benchmark_parcel_merging_LDADD = 			\
	$(top_builddir)/src/utils/libepic_utils.la 	\
	$(top_builddir)/src/mpi/libepic_mpi.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

clean-local:
	rm -f *.mod *.smod
