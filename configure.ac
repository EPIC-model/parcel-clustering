AC_INIT([parcel-clustering],
        [1.0.0],
        [mf248@st-andrews.ac.uk],
        [],
        [https://github.com/EPIC-model/parcel-clustering])
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_PREREQ([2.69])

# reset default options:
# (see also https://www.gnu.org/savannah-checkouts/gnu/autoconf/manual/autoconf-2.72/autoconf.pdf)
: ${CFLAGS=""}
: ${FCFLAGS=""}
: ${CXXFLAGS=""}

AC_PROG_FC([ftn mpifort])
AC_PROG_CC([cc mpicc])
AC_LANG(Fortran)
AM_PROG_AR
LT_INIT

# change file extension from *.f to *.f90
# (important for library tests since it autogenerates a file conftest.f90)
AC_FC_SRCEXT(f90)
AC_FC_PP_SRCEXT(f90)

AC_CONFIG_MACRO_DIRS([m4])
AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([
    Makefile
    cpp/Makefile
    src/Makefile
    src/utils/Makefile
    pytools/Makefile
])

# 30 May 2023
# https://www.gnu.org/software/autoconf-archive/ax_compiler_vendor.html
# https://www.gnu.org/software/autoconf-archive/ax_check_compile_flag.html
AX_COMPILER_VENDOR

if test "$ax_cv_fc_compiler_vendor" = "unknown"; then
    AC_MSG_ERROR([Could not deduce compiler vendor.])
fi

# gfortran compiler flags
if test "$ax_cv_fc_compiler_vendor" = "gnu"; then
    AX_CHECK_COMPILE_FLAG([-cpp], [FCFLAGS="$FCFLAGS -cpp"])
    AX_CHECK_COMPILE_FLAG([-mcmodel=large], [FCFLAGS="$FCFLAGS -mcmodel=large"])
    AX_CHECK_COMPILE_FLAG([-std=f2018], [FCFLAGS="$FCFLAGS -std=f2018"])
    AX_CHECK_COMPILE_FLAG([-fall-intrinsics], [FCFLAGS="$FCFLAGS -fall-intrinsics"])
    AX_CHECK_COMPILE_FLAG([-fdefault-real-8], [FCFLAGS="$FCFLAGS -fdefault-real-8"])
    AX_CHECK_COMPILE_FLAG([-fdefault-double-8], [FCFLAGS="$FCFLAGS -fdefault-double-8"])
    AX_CHECK_COMPILE_FLAG([-Wno-unused-dummy-argument], [FCFLAGS="$FCFLAGS -Wno-unused-dummy-argument"])
    AX_CHECK_COMPILE_FLAG([-fcray-pointer], [FCFLAGS="$FCFLAGS -fcray-pointer"])
fi

# intel compiler flags
if test "$ax_cv_fc_compiler_vendor" = "intel"; then
    AX_CHECK_COMPILE_FLAG([-std18], [FCFLAGS="$FCFLAGS -std18"])
    AX_CHECK_COMPILE_FLAG([-fp-model=source], [FCFLAGS="$FCFLAGS -fp-model=source"])
    AX_CHECK_COMPILE_FLAG([-real-size 8], [FCFLAGS="$FCFLAGS -real-size 8"])
    AX_CHECK_COMPILE_FLAG([-double-size 64], [FCFLAGS="$FCFLAGS -double-size 64"])
    # Suppress harmless warnings:
    # > warning #8100: The actual argument is an array section or assumed-shape array,
    # > corresponding dummy argument that has either the VOLATILE or ASYNCHRONOUS attribute shall be an assumed-shape array.   [SEND_BUF]
    # >       call MPI_Isend(send_buf(1:send_size)
    AX_CHECK_COMPILE_FLAG([-diag-disable 8100], [FCFLAGS="$FCFLAGS -diag-disable 8100"])
fi

# Cray Fortran compiler:
if test "$ax_cv_fc_compiler_vendor" = "cray"; then
    AX_CHECK_COMPILE_FLAG([-s integer32], [FCFLAGS="$FCFLAGS -s integer32"])
    AX_CHECK_COMPILE_FLAG([-s real64], [FCFLAGS="$FCFLAGS -s real64"])
    # Suppress harmless warnings:
    #
    # The following text is output from 'explain ftn-878':
    # > Warning : A module named "%s" has already been directly or indirectly use
    # > associated into this scope.
    # >
    # > A module name is a global name and by definition, all global names must be
    # > unique.  The compiler detected more than one module with the same name in
    # > this scope and assumes that they are the same module.  This can happen
    # > because a module can be directly or indirectly use associated into a scoping
    # > unit.  An example would be:
    # >
    # >    MODULE ABC
    # >    INTEGER I
    # >    END MODULE
    # >
    # >    MODULE CDE
    # >    USE ABC
    # >    END MODULE
    # >
    # >    PROGRAM BBB
    # >    USE CDE
    # >    USE ABC  ! ABC is already in this scope via CDE
    # >    END PROGRAM
    #
    # The following text is output from 'explain ftn-1077':
    # > Warning : This compilation contains OpenMP directives.  -h omp is not active
    # > so the directives are ignored.
    # >
    # > OpenMP directives are only recognized when -h omp is specified.  This is not
    # > the default mode.   The compiler has found one or more OpenMP directives.
    # > These directives have been ignored.  This message is only issued once per
    # > compilation.  This message can be disabled with -M1077.  Note that -h
    # > profile_generate overrides -h omp, so if both are specified, -h omp will be
    # > turned off.
	AX_CHECK_COMPILE_FLAG([-M 878,1077], [FCFLAGS="$FCFLAGS -M 878,1077"])
fi

#######################################################################################
##
## "--with" flags
##


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NetCDF C

## 11 March 2021
## https://www.gnu.org/software/autoconf/manual/autoconf-2.60/html_node/External-Software.html
AC_ARG_WITH([netcdf-c],
            [AS_HELP_STRING([--with-netcdf-c=<path>], [path to NetCDF C library @<:@default=$NETCDF_C_DIR@:>@])],
            [NETCDF_C_DIR=$withval],
            [])

if test -n "${NETCDF_C_DIR}"; then
    CFLAGS="$CFLAGS -I$NETCDF_C_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_C_DIR/lib"
else
    AC_MSG_ERROR([No NETCDF_C_DIR environment variable.])
fi

AC_LANG_PUSH([C])
AC_SEARCH_LIBS([nc_create], [netcdf], [], [], [])
AC_MSG_CHECKING([whether we can compile a NetCDF C program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([#include <netcdf.h>], [
    int ncid;
    nc_create("dummy.nc", NC_CLOBBER, &ncid);
    ])],
[netcdf_c_found=yes],
[netcdf_c_found=no])
AC_MSG_RESULT([$netcdf_c_found])
AC_LANG_POP([C])

if test "x$netcdf_c_found" = "xno"; then
    AC_MSG_ERROR([Cannot compile a NetCDF C program])
fi


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# NetCDF Fortran

AC_ARG_WITH([netcdf-fortran],
            [AS_HELP_STRING([--with-netcdf-fortran=<path>],
            [path to NetCDF Fortran library @<:@default=$NETCDF_FORTRAN_DIR@:>@])],
            [NETCDF_FORTRAN_DIR=$withval],
            [])

if test -n "${NETCDF_FORTRAN_DIR}"; then
    FCFLAGS="$FCFLAGS -I$NETCDF_FORTRAN_DIR/include"
    LDFLAGS="$LDFLAGS -L$NETCDF_FORTRAN_DIR/lib"
else
    AC_MSG_ERROR([No NETCDF_FORTRAN_DIR environment variable.])
fi


# 15 April 2023
#http://git.savannah.gnu.org/gitweb/?p=autoconf-archive.git;a=blob_plain;f=m4/ax_prog_fc_mpi.m4
# We cannot use AC_SEARCH_LIBS
for lib in NONE netcdff; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([for function nf90_open])
    else
    AC_MSG_CHECKING([for function nf90_open in -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [
           use netcdf
            implicit none
            integer :: ncerr, ncid
            ncerr = nf90_open("dummy.nc", NF90_NOWRITE, ncid)])],
        [netcdf_found=yes],
        [netcdf_found=no])
    AC_MSG_RESULT($netcdf_found)
    if test "x$netcdf_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done

AC_MSG_CHECKING([whether we can compile a NetCDF Fortran program])
AC_LINK_IFELSE(
[AC_LANG_PROGRAM([], [
    use netcdf
    implicit none])],
[netcdf_fortran_found=yes],
[netcdf_fortran_found=no])
AC_MSG_RESULT([$netcdf_fortran_found])

if test "x$netcdf_fortran_found" = "xno"; then
    AC_MSG_ERROR([Cannot compile a NetCDF Fortran program])
fi


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# MPI

AC_ARG_WITH([mpi],
            [AS_HELP_STRING([--with-mpi=<path>], [path to mpi library @<:@default=$MPI_DIR@:>@])],
            [MPI_DIR=$withval],
            [])

if test -n "${MPI_DIR}"; then
    # mod file may also be in $MPI_DIR/lib
    AC_MSG_NOTICE([Found MPI_DIR environment. Appending include and lib path to FCFLAGS and LDFLAGS.])
    FCFLAGS="$FCFLAGS -I$MPI_DIR/include -I$MPI_DIR/lib"
    LDFLAGS="$LDFLAGS -L$MPI_DIR/lib"
fi

AC_SEARCH_LIBS([MPI_Init],
               [mpi mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90 mpifort_cray fmpich mpifort],
               [],
               [])

AC_SEARCH_LIBS([MPI_Win_lock],
               [mpi mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90 mpifort_cray fmpich mpifort],
               [],
               [])

AC_SEARCH_LIBS([MPI_Win_flush],
               [mpi mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90 mpifort_cray fmpich mpifort],
               [],
               [])

AC_SEARCH_LIBS([MPI_Wtime],
               [mpi mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90 mpifort_cray fmpich mpifort],
               [],
               [])

if test "$ax_cv_fc_compiler_vendor" = "gnu"; then
for lib in NONE mpi mpi_usempif08 mpi_usempi_ignore_tkr mpi_mpifh mpichf90 mpifort_cray fmpich mpifort; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([for module mpi_f08])
    else
    AC_MSG_CHECKING([for module mpi_f08 in -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
        [mpi_found=yes],
        [mpi_found=no])
    AC_MSG_RESULT([$mpi_found])

    if test "x$mpi_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done
fi

AC_MSG_CHECKING([for MPI Fortran library])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
    [mpi_found=yes],
    [mpi_found=no])
AC_MSG_RESULT([$mpi_found])

if test "x$mpi_found" = "xno"; then
    AC_MSG_ERROR([Cannot find MPI Fortran library])
fi

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check if MPI_INTEGER8 exists
AC_MSG_CHECKING([for MPI_INTEGER8])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
        type(MPI_Datatype) :: dtype
        call MPI_Init(ierr)
        dtype = MPI_INTEGER8
        call MPI_Finalize(ierr)])],
    [mpi_integer8_found=yes],
    [mpi_integer8_found=no])
AC_MSG_RESULT([$mpi_integer8_found])

if test "x$mpi_integer8_found" = "xyes"; then
    AX_CHECK_COMPILE_FLAG([-DENABLE_MPI_INTEGER8], [FCFLAGS="$FCFLAGS -DENABLE_MPI_INTEGER8"])
fi

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Check if we can write e.g. "type(MPI_Comm) :: comm = MPI_COMM_NULL"
# E.g. with ifort we get the following error:
# This symbol must be a defined parameter, an enumerator, or an argument of an inquiry function that evaluates to a compile-time constant.
AC_MSG_CHECKING([if we can compile with direct assignment of MPI_COMM_NULL and MPI_OP_NULL ])
AC_LINK_IFELSE(
    [AC_LANG_PROGRAM([], [
        use mpi_f08
        implicit none
        integer :: ierr
	type(MPI_Comm) :: comm = MPI_COMM_NULL
	type(MPI_Op) :: op = MPI_OP_NULL
        call MPI_Init(ierr)
        call MPI_Finalize(ierr)])],
    [null_assign_works=yes],
    [null_assign_works=no])
AC_MSG_RESULT([$null_assign_works])

if test "x$null_assign_works" = "xyes"; then
    AX_CHECK_COMPILE_FLAG([-DNULL_ASSIGNMENT_WORKS], [FCFLAGS="$FCFLAGS -DNULL_ASSIGNMENT_WORKS"])
fi

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# OpenSHMEM

AC_LANG_PUSH([C])

AC_SEARCH_LIBS([shmem_init],
            [oshmem],
            [],
            [])

AC_SEARCH_LIBS([shmem_finalize],
            [oshmem],
            [],
            [])

AC_SEARCH_LIBS([shmem_putmem],
            [oshmem],
            [],
            [])

AC_SEARCH_LIBS([shmem_getmem],
            [oshmem],
            [],
            [])


AC_SEARCH_LIBS([shmem_put32],
            [oshmem],
            [],
            [])

AC_SEARCH_LIBS([shmem_get32],
            [oshmem],
            [],
            [])

for lib in NONE oshmem; do
    save_LIBS=$LIBS
    if test x"$lib" = xNONE; then
    AC_MSG_CHECKING([whether we can compile a SHMEM C program])
    else
    AC_MSG_CHECKING([whether we can compile a SHMEM C program using -l$lib])
    LIBS="-l$lib $LIBS"
    fi
    AC_LINK_IFELSE(
        [AC_LANG_PROGRAM([#include <shmem.h>], [
        shmem_init();
        int rank = shmem_my_pe();
        int nranks = shmem_n_pes();
        shmem_finalize();])],
        [openshmem_found=yes],
        [openshmem_found=no])
    AC_MSG_RESULT([$openshmem_found])

    if test "x$openshmem_found" = "xyes"; then
    break;
    fi
    LIBS=$save_LIBS
done

AC_LANG_POP([C])

AC_MSG_CHECKING([whether we are enabling SHMEM])
if test "x$openshmem_found" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_SHMEM"
else
    AC_MSG_RESULT([no])
fi

AM_CONDITIONAL([ENABLE_SHMEM], [test "$openshmem_found" = "yes"])

if test "x$openshmem_found" = "xyes"; then
    AC_MSG_CHECKING([size of Fortran logical])
    # (accessed 5 March 2025)
    # https://www.gnu.org/software/autoconf/manual/autoconf-2.63/html_node/Runtime.html
    AC_RUN_IFELSE([
        AC_LANG_PROGRAM([],
                        [use, intrinsic :: iso_c_binding
                        implicit none
                        logical :: var
                        open(unit=12, file='conftest.out', status='replace')
                        write(12, fmt='(I2)') 8*c_sizeof(var) ! write size in bits
                        close(12)])],
        [logical_size=$(cat conftest.out)],
        [AC_MSG_ERROR([Could not figure out size of Fortran logical.])])

    AC_MSG_RESULT([$logical_size bit])
    if test $logical_size == 32; then
        FCFLAGS="$FCFLAGS -DENABLE_SHMEM_PUT_GET_32"
    fi
fi

#######################################################################################
##
## "--enable" flags
##

ENABLE_COARRAY='no'
AC_ARG_ENABLE([coarray],
              [AS_HELP_STRING([--enable-coarray], [enable coarray mode (default=no)])],
              [ENABLE_COARRAY=$enableval])

AM_CONDITIONAL([ENABLE_COARRAY], [test "$ENABLE_COARRAY" = "yes"])

AC_MSG_CHECKING([whether we are enabling Fortran Coarray])
if test "x$ENABLE_COARRAY" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_COARRAY"
    AX_CHECK_COMPILE_FLAG([-hcaf], [FCFLAGS="$FCFLAGS -hcaf"])
else
    AC_MSG_RESULT([no])
fi


##
## Python
##
## # 23 March 2021
## https://askubuntu.com/questions/29370/how-to-check-if-a-command-succeeded
## https://gitlab.psi.ch/H5hut/src/-/blob/master/configure.ac
##

ENABLE_PYTHON_CHECK='no'
AC_ARG_ENABLE([python],
              [AS_HELP_STRING([--enable-python-check], [enable Python check (default=no)])],
              [ENABLE_PYTHON_CHECK=$enableval])

AM_CONDITIONAL([ENABLE_PYTHON_CHECK], [test "$ENABLE_PYTHON_CHECK" = "yes"])

if test "x$ENABLE_PYTHON_CHECK" = "xyes"; then
    AM_PATH_PYTHON([3.11])
    #
    # argparse
    #
    AC_MSG_CHECKING([for python argparse library])
    `${PYTHON} -c 'import argparse' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install argparse with "conda install conda-forge::argparse".])
        exit 1
    fi

    #
    # numpy
    #
    AC_MSG_CHECKING([for python numpy library])
    `${PYTHON} -c 'import numpy' 2> config.log`
    if test "X$?" = "X0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install numpy with "conda install conda-forge::numpy".])
        exit 1
    fi

    #
    # pandas
    #
    AC_MSG_CHECKING([for python pandas library])
    `${PYTHON} -c 'import pandas' 2> config.log`
    if test "x$?" = "x0" ; then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_NOTICE([Install pandas with "conda install conda-forge::pandas".])
        exit 1
    fi
fi

ENABLE_VERBOSE='no'
AC_ARG_ENABLE([verbose],
              [AS_HELP_STRING([--enable-verbose], [enable verbose mode (default=no)])],
              [ENABLE_VERBOSE=$enableval])

AM_CONDITIONAL([ENABLE_VERBOSE], [test "$ENABLE_VERBOSE" = "yes"])

AC_MSG_CHECKING([whether we are compiling in verbose mode])
if test "x$ENABLE_VERBOSE" = "xyes"; then
    AC_MSG_RESULT([yes])
    FCFLAGS="$FCFLAGS -DENABLE_VERBOSE"
else
    AC_MSG_RESULT([no])
fi


ENABLE_DEBUG='no'
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug], [compile in debug mode (default=no)])],
              [ENABLE_DEBUG=$enableval])

AM_CONDITIONAL([ENABLE_DEBUG], [test "$ENABLE_DEBUG" = "yes"])

AC_MSG_CHECKING([whether we are compiling in debug mode])
if test "x$ENABLE_DEBUG" = "xyes"; then
    AC_MSG_RESULT([yes])
    AX_CHECK_COMPILE_FLAG([-g], [FCFLAGS="$FCFLAGS -g"])
    AX_CHECK_COMPILE_FLAG([-O0], [FCFLAGS="$FCFLAGS -O0"])

    # gfortran compiler flags
    if test "$ax_cv_fc_compiler_vendor" = "gnu"; then
        AX_CHECK_COMPILE_FLAG([-Wall], [FCFLAGS="$FCFLAGS -Wall"])
        AX_CHECK_COMPILE_FLAG([-Wuninitialized], [FCFLAGS="$FCFLAGS -Wuninitialized"])
        AX_CHECK_COMPILE_FLAG([-Wmaybe-uninitialized], [FCFLAGS="$FCFLAGS -Wmaybe-uninitialized"])
        AX_CHECK_COMPILE_FLAG([-Werror], [FCFLAGS="$FCFLAGS -Werror"])
        AX_CHECK_COMPILE_FLAG([-fcheck=all], [FCFLAGS="$FCFLAGS -fcheck=all"])
        AX_CHECK_COMPILE_FLAG([-fbounds-check], [FCFLAGS="$FCFLAGS -fbounds-check"])
        AX_CHECK_COMPILE_FLAG([-fbacktrace], [FCFLAGS="$FCFLAGS -fbacktrace"])
        AX_CHECK_COMPILE_FLAG([-ffpe-trap=denormal,invalid,zero,overflow,underflow],
                              [FCFLAGS="$FCFLAGS -ffpe-trap=denormal,invalid,zero,overflow,underflow"])
    fi

    # intel compiler flags
    if test "$ax_cv_fc_compiler_vendor" = "intel"; then
        AX_CHECK_COMPILE_FLAG([-warn all], [FCFLAGS="$FCFLAGS -Wall"])
        AX_CHECK_COMPILE_FLAG([-warn error], [FCFLAGS="$FCFLAGS -warn error"])
        AX_CHECK_COMPILE_FLAG([-debug full], [FCFLAGS="$FCFLAGS -debug full"])
    fi
else
    AC_MSG_RESULT([no])
    AX_CHECK_COMPILE_FLAG([-DNDEBUG], [FCFLAGS="$FCFLAGS -DNDEBUG"])
    AX_CHECK_COMPILE_FLAG([-O2],  [FCFLAGS="$FCFLAGS -O2"])

    # gfortran compiler flags
    if test "$ax_cv_fc_compiler_vendor" = "gnu"; then
        AX_CHECK_COMPILE_FLAG([-funroll-all-loops], [FCFLAGS="$FCFLAGS -funroll-all-loops"])
    fi

    # intel compiler flags
    if test "$ax_cv_fc_compiler_vendor" = "intel"; then
        AX_CHECK_COMPILE_FLAG([-funroll-loops], [FCFLAGS="$FCFLAGS -funroll-loops"])
        AX_CHECK_COMPILE_FLAG([-ipo], [FCFLAGS="$FCFLAGS -ipo"])
    fi

    # Cray fortran compiler flags
    if test "$ax_cv_fc_compiler_vendor" = "cray"; then
        AX_CHECK_COMPILE_FLAG([-h ipa2], [FCFLAGS="$FCFLAGS -h ipa2"])
        AX_CHECK_COMPILE_FLAG([-h fp2], [FCFLAGS="$FCFLAGS -h fp2"])
        AX_CHECK_COMPILE_FLAG([-h unroll2],  [FCFLAGS="$FCFLAGS -h unroll2"])
    fi
fi


ENABLE_OPENMP='no'
AC_ARG_ENABLE([openmp],
              [AS_HELP_STRING([--enable-openmp], [enable OpenMP (default=no)])],
              [ENABLE_OPENMP=$enableval])

AM_CONDITIONAL([ENABLE_OPENMP], [test "$ENABLE_OPENMP" = "yes"])

AC_MSG_CHECKING([whether we are enabling OpenMP])
if test "x$ENABLE_OPENMP" = "xyes"; then
    AC_MSG_RESULT([yes])
    # gfortran compiler flags
    AX_CHECK_COMPILE_FLAG([-fopenmp], [FCFLAGS="$FCFLAGS -fopenmp -DENABLE_OPENMP"])

    # intel compiler flags
    AX_CHECK_COMPILE_FLAG([-qopenmp], [FCFLAGS="$FCFLAGS -qopenmp -DENABLE_OPENMP"])
else
    AC_MSG_RESULT([no])
fi

AC_PROG_CXX([CC mpic++])
AC_LANG_PUSH([C++])
AX_CHECK_COMPILE_FLAG([-std=c++17], [CXXFLAGS="$CXXFLAGS -std=c++17"])
AC_LANG_POP([C++])


AC_OUTPUT

echo ""
echo "Compilers:"
echo "----------"
echo "Fortran compiler: $(which $FC)"
echo "C compiler:       $(which $CC)"
echo "C++ compiler:     $(which $CXX)"
echo ""
echo "Compiler flags:"
echo "---------------"
echo "FCLAGS: $FCFLAGS"
echo ""
echo "CFLAGS: $CFLAGS"
echo ""
echo "CXXFLAGS: $CXXFLAGS"
echo ""
