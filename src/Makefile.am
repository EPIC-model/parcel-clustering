SUBDIRS = utils

# this line is needed such that config.h is found
AM_FCFLAGS = -I $(top_builddir)/src/		\
	     -I $(top_builddir)/src/utils/

bin_PROGRAMS = benchmark_random			\
               benchmark_read

if ENABLE_PYTHON
bin_PROGRAMS += benchmark_verify
endif

benchmark_random_SOURCES = 			\
	mpi/mpi_tags.f90                        \
        mpi/mpi_datatypes.f90                   \
        mpi/mpi_ops.f90                         \
        mpi/mpi_environment.f90                 \
        mpi/mpi_utils.f90                       \
        mpi/mpi_layout.f90                      \
        mpi/mpi_collectives.f90                 \
        mpi/mpi_timer.f90                       \
        config.f90                              \
        netcdf/netcdf_utils.f90                 \
        netcdf/netcdf_writer.f90                \
        netcdf/netcdf_reader.f90                \
        netcdf/netcdf_timings.f90		\
        options.f90                             \
        parameters.f90                          \
        armanip.f90                             \
        parcels/parcel_ops.f90                  \
        parcels/parcel_container.f90            \
        parcels/parcel_ellipsoid.f90            \
        parcels/parcels.f90                     \
        parcels/parcel_bc.f90                   \
        parcels/parcel_mpi.f90                  \
        parcels/parcel_split.f90                \
        parcels/parcel_init.f90                 \
        parcels/parcel_nearest_graph.f90        \
        parcels/parcel_nearest_p2p_graph.f90    \
        parcels/parcel_nearest_rma_graph.f90    \
        parcels/parcel_nearest_shmem_graph.f90  \
        parcels/parcel_nearest.f90              \
        parcels/parcel_merge.f90                \
        parcels/parcel_netcdf.f90               \
        utils.f90                               \
        serial/parcel_nearest_serial.f90        \
        serial/parcel_merge_serial.f90          \
	benchmark_random.f90

benchmark_random_LDADD = $(top_builddir)/src/utils/libcluster_utils.la

benchmark_read_SOURCES = 			\
	mpi/mpi_tags.f90                        \
        mpi/mpi_datatypes.f90                   \
        mpi/mpi_ops.f90                         \
        mpi/mpi_environment.f90                 \
        mpi/mpi_utils.f90                       \
        mpi/mpi_layout.f90                      \
        mpi/mpi_collectives.f90                 \
        mpi/mpi_timer.f90                       \
        config.f90                              \
        netcdf/netcdf_utils.f90                 \
        netcdf/netcdf_writer.f90                \
        netcdf/netcdf_reader.f90                \
        netcdf/netcdf_timings.f90		\
        options.f90                             \
        parameters.f90                          \
        armanip.f90                             \
        parcels/parcel_ops.f90                  \
        parcels/parcel_container.f90            \
        parcels/parcel_ellipsoid.f90            \
        parcels/parcels.f90                     \
        parcels/parcel_bc.f90                   \
        parcels/parcel_mpi.f90                  \
        parcels/parcel_split.f90                \
        parcels/parcel_init.f90                 \
        parcels/parcel_nearest_graph.f90        \
        parcels/parcel_nearest_p2p_graph.f90    \
        parcels/parcel_nearest_rma_graph.f90    \
        parcels/parcel_nearest_shmem_graph.f90  \
        parcels/parcel_nearest.f90              \
        parcels/parcel_merge.f90                \
        parcels/parcel_netcdf.f90               \
        utils.f90                               \
        serial/parcel_nearest_serial.f90        \
        serial/parcel_merge_serial.f90          \
	benchmark_read.f90

benchmark_read_LDADD = $(top_builddir)/src/utils/libcluster_utils.la

if ENABLE_PYTHON
benchmark_verify_SOURCES = 			\
	mpi/mpi_tags.f90                        \
        mpi/mpi_datatypes.f90                   \
        mpi/mpi_ops.f90                         \
        mpi/mpi_environment.f90                 \
        mpi/mpi_utils.f90                       \
        mpi/mpi_layout.f90                      \
        mpi/mpi_collectives.f90                 \
        mpi/mpi_timer.f90                       \
        config.f90                              \
        netcdf/netcdf_utils.f90                 \
        netcdf/netcdf_writer.f90                \
        netcdf/netcdf_reader.f90                \
        netcdf/netcdf_timings.f90		\
        options.f90                       	\
        parameters.f90                    	\
        armanip.f90                       	\
        parcels/parcel_ops.f90                  \
        parcels/parcel_container.f90            \
        parcels/parcel_ellipsoid.f90            \
        parcels/parcels.f90                     \
        parcels/parcel_bc.f90                   \
        parcels/parcel_mpi.f90                  \
        parcels/parcel_split.f90                \
        parcels/parcel_init.f90                 \
        parcels/parcel_nearest_graph.f90        \
        parcels/parcel_nearest_p2p_graph.f90    \
        parcels/parcel_nearest_rma_graph.f90    \
        parcels/parcel_nearest_shmem_graph.f90  \
        parcels/parcel_nearest.f90              \
        parcels/parcel_merge.f90                \
        parcels/parcel_netcdf.f90               \
        utils.f90                               \
        serial/parcel_nearest_serial.f90        \
        serial/parcel_merge_serial.f90		\
	benchmark_verify.f90

benchmark_verify_LDADD = $(top_builddir)/src/utils/libcluster_utils.la
endif

clean-local:
	rm -f *.mod *.smod
