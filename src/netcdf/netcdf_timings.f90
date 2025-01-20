module netcdf_timings
    use mpi_timer, only : timings
    use mpi_collectives, only : mpi_blocking_reduce
    use config, only : package_version, cf_version
    use mpi_utils, only : mpi_exit_on_error
    use netcdf_utils
    use netcdf_writer
    use netcdf_reader, only : get_num_steps &
                            , get_var_id
    use iomanip, only : int2string
    implicit none

    private

    integer :: ncid = -1
    integer :: n_writes = 0
    integer :: unlimdimid = -1


    type(netcdf_info), allocatable :: nc_dset(:)

    public :: write_netcdf_timings

contains

    subroutine create_file(ncfname)
        character(len=*), intent(in) :: ncfname
        integer                      :: n_timers, n, m

        call create_netcdf_file(ncfname,            &
                                overwrite=.false.,  &
                                ncid=ncid,          &
                                l_serial=.true.)

        call write_netcdf_info(ncid=ncid,                    &
                               version_tag=package_version,  &
                               file_type='timing',           &
                               cf_version=cf_version)

        call define_netcdf_dimension(ncid=ncid,                 &
                                     name='t',                  &
                                     dimsize=NF90_UNLIMITED,    &
                                     dimid=unlimdimid)

        n_timers = size(timings)
        allocate(nc_dset(2 * n_timers))

        ! max wall time
        do n = 1, n_timers
            nc_dset(n)%name = "wtime-" // trim(int2string(n))
            nc_dset(n)%dtype = NF90_DOUBLE
            call define_netcdf_dataset(ncid=ncid,                       &
                                       name=nc_dset(n)%name,            &
                                       long_name=trim(timings(n)%name), &
                                       std_name="",                     &
                                       unit="s",                        &
                                       dtype=nc_dset(n)%dtype,          &
                                       dimids=(/unlimdimid/),           &
                                       varid=nc_dset(n)%varid)
        enddo

        ! max number of calls
        do m = 1, n_timers
            n = n_timers + m
            nc_dset(n)%name = "ncalls-" // trim(int2string(m))
            nc_dset(n)%dtype = NF90_INT
            call define_netcdf_dataset(ncid=ncid,                       &
                                       name=nc_dset(n)%name,            &
                                       long_name=trim(timings(m)%name), &
                                       std_name="",                     &
                                       unit="1",                        &
                                       dtype=nc_dset(n)%dtype,          &
                                       dimids=(/unlimdimid/),           &
                                       varid=nc_dset(n)%varid)
        enddo

        call close_definition(ncid)

        call close_netcdf_file(ncid, l_serial=.true.)

    end subroutine create_file

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine read_file_content(ncfname)
        character(len=*), intent(in) :: ncfname
        integer                      :: n, m, n_vars, n_timers

        call open_netcdf_file(ncfname, NF90_NOWRITE, ncid, l_serial=.true.)
        call get_num_steps(ncid, n_writes)

        ncerr = nf90_inquire(ncid, nVariables=n_vars, unlimitedDimID=unlimdimid)

        call check_netcdf_error("Unable to retrieve netCDF info")

        n_timers = size(timings)

        if (n_vars /= 2 * n_timers) then
            call mpi_exit_on_error("Number of timers disagree.")
        endif

        allocate(nc_dset(n_vars))

        do n = 1, n_timers
            nc_dset(n)%name = "wtime-" // trim(int2string(n))
            call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
        enddo

        do m = 1, n_timers
            n = n_timers + m
            nc_dset(n)%name = "ncalls-" // trim(int2string(m))
            call get_var_id(ncid, nc_dset(n)%name, nc_dset(n)%varid)
        enddo

        call close_netcdf_file(ncid, l_serial=.true.)

    end subroutine read_file_content

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine write_netcdf_timings(ncfname)
        character(len=*), intent(in)  :: ncfname
        logical                       :: l_exist
        integer                       :: n, m, n_timers
        double precision, allocatable :: wbuffer(:)
        integer, allocatable          :: nbuffer(:)

        n_timers = size(timings)

        allocate(wbuffer(1:n_timers))
        allocate(nbuffer(1:n_timers))
        wbuffer(1:n_timers) = timings(1:n_timers)%wall_time
        nbuffer(1:n_timers) = timings(1:n_timers)%n_calls

        call mpi_blocking_reduce(wbuffer, MPI_MAX, world)
        call mpi_blocking_reduce(nbuffer, MPI_MAX, world)

        if (world%rank /= world%root) then
            deallocate(wbuffer)
            deallocate(nbuffer)
            return
        endif

        call exist_netcdf_file(ncfname, l_exist)

        if (.not. l_exist) then
            call create_file(ncfname)
        else
            call read_file_content(ncfname)
        endif

        call open_netcdf_file(ncfname, NF90_WRITE, ncid, l_serial=.true.)

        n_writes = n_writes + 1

        do n = 1, n_timers
            m = find_loc("wtime-", n)

            if (m == 0) then
                print *, "NetCDF variable wtime-" // trim(int2string(n)) // " not found."
                cycle
            endif
#ifndef NDEBUG
            if (nc_dset(m)%name /= "wtime-" // trim(int2string(n))) then
                call mpi_exit_on_error("Timer names disagree.")
            endif
#endif
            call write_netcdf_scalar(ncid, nc_dset(m)%varid, wbuffer(n), &
                                     n_writes, l_serial=.true.)
        enddo

        do n = 1, n_timers
            m = find_loc("ncalls-", n)

            if (m == 0) then
                print *, "NetCDF variable ncalls-" // trim(int2string(n)) // " not found."
                cycle
            endif

#ifndef NDEBUG
            if (nc_dset(m)%name /= "ncalls-" // trim(int2string(n))) then
                call mpi_exit_on_error("Timer names disagree.")
            endif
#endif
            call write_netcdf_scalar(ncid, nc_dset(m)%varid, nbuffer(n), &
                                     n_writes, l_serial=.true.)
        enddo

        call close_netcdf_file(ncid, l_serial=.true.)

        deallocate(wbuffer)
        deallocate(nbuffer)

    end subroutine write_netcdf_timings

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Replaces the intrinsic Fortran routine which does not behave with CCE like the GCC version
    ! m = findloc(nc_dset(1:n_timers)%name, value = "ncalls-" // trim(int2string(n)), dim=1)
    function find_loc(s, n) result(r)
        character(*), intent(in) :: s
        integer,      intent(in) :: n
        integer                  :: m, r

        r = 0
        do m = 1, size(nc_dset(:)%name)
            if (nc_dset(m)%name == s // trim(int2string(n))) then
                r = m
                exit
            endif
        enddo
    end function find_loc


end module netcdf_timings
