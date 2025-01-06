AM_FCFLAGS =				\
	-I $(top_builddir)/src/mpi	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/netcdf

lib_LTLIBRARIES = libparcels.la
libparcels_la_SOURCES =			\
	parcel_ops.f90			\
	parcel_container.f90 		\
	parcel_ellipsoid.f90		\
	parcels.f90			\
	parcel_bc.f90                   \
	parcel_mpi.f90			\
	parcel_split.f90		\
	parcel_init.f90			\
	parcel_nearest_graph.f90	\
	parcel_nearest_p2p_graph.f90	\
	parcel_nearest_rma_graph.f90	\
	parcel_nearest_shmem_graph.f90	\
	parcel_nearest.f90		\
	parcel_merge.f90		\
	parcel_netcdf.f90
