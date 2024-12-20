# this line is needed such that config.h is found
AM_FCFLAGS = -I $(top_builddir)/src/		\
	     -I $(top_builddir)/src/mpi		\
	     -I $(top_builddir)/src/netcdf

lib_LTLIBRARIES = libepic_utils.la
libepic_utils_la_SOURCES = 	\
	dimensions.f90		\
	datatypes.f90		\
	armanip.f90		\
	iomanip.f90		\
	constants.f90		\
	options.f90		\
	parameters.f90		\
	linalg.f90		\
	scherzinger.f90

