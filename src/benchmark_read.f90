program benchmark_read
    use constants, only : pi, zero, one, f12, f23, twopi
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_init, only : parcel_default
    use parcel_nearest, only : tree
    use parcel_mpi, only : parcel_communicate
    use mpi_environment
    use mpi_layout
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_ops, only : MPI_SUM_64BIT
    use mpi_utils, only : mpi_stop
#ifdef ENABLE_COARRAY
    use mpi_utils, only : mpi_print
#endif
    use utils, only : register_all_timers
    use parcel_merging
    use parcel_netcdf
    use netcdf_utils
    use netcdf_reader
    use iomanip, only : zfill
#ifndef ENABLE_COARRAY
    use parcel_nearest_p2p_graph, only : p2p_graph_t
    use parcel_nearest_rma_graph, only : rma_graph_t
#ifdef ENABLE_SHMEM
    use parcel_nearest_shmem_graph, only : shmem_graph_t
#endif
#endif
    use netcdf_timings
    implicit none

    character(64)    :: basename
    character(512)   :: fname, ncfname
    integer          :: ncid, n, m, niter
    integer          :: ncells(3), offset, nfiles
    character(len=5) :: comm_type ! p2p, rma, shmem or caf
    logical          :: l_subcomm

    call mpi_env_initialise

    call register_all_timers

    parcel%lambda_max = 4.0d0
    parcel%min_vratio = 20.0d0
    parcel%size_factor = 1.0d0

    call parse_command_line

    fname = trim(basename) // '_' // zfill(offset) // '_parcels.nc'
    call open_netcdf_file(trim(fname), NF90_NOWRITE, ncid)

    call get_netcdf_box(ncid, lower, extent, ncells)
    call close_netcdf_file(ncid)

    nx = ncells(1)
    ny = ncells(2)
    nz = ncells(3)

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call parcels%allocate(max_num_parcels)

#ifndef ENABLE_COARRAY
    select case(comm_type)
        case ('p2p')
            allocate(p2p_graph_t :: tree)
        case ('rma')
            allocate(rma_graph_t :: tree)
#ifdef ENABLE_SHMEM
        case ('shmem')
            allocate(shmem_graph_t :: tree)
#endif
        case default
            call mpi_stop("Communication layer not available!")
    end select
#endif

    call tree%initialise(max_num_parcels, l_subcomm)

    call tree%register_timer

    do n = 0, niter - 1

        m = mod(n, nfiles) + offset

        fname = trim(basename) // '_' // zfill(m) // '_parcels.nc'
        if (world%rank == world%root) then
           print *, "Read:", trim(fname)
        endif
        call read_netcdf_parcels(fname)

        parcels%total_num = 0

        call MPI_Allreduce(parcels%local_num, &
                           parcels%total_num, &
                           1,                 &
                           MPI_INTEGER_64BIT, &
                           MPI_SUM_64BIT,     &
                           world%comm,        &
                           world%err)

        if (world%rank == world%root) then
            print *, "Number of parcels before merging:", parcels%total_num
        endif

        call parcel_merge

        parcels%total_num = 0
        call MPI_Allreduce(parcels%local_num, &
                           parcels%total_num, &
                           1,                 &
                           MPI_INTEGER_64BIT, &
                           MPI_SUM_64BIT,     &
                           world%comm,        &
                           world%err)


        if (world%rank == world%root) then
            print *, "Number of parcels after merging:", parcels%total_num
        endif
    enddo

    call parcels%deallocate

    call write_netcdf_timings(trim(ncfname))

    call tree%finalise

    call mpi_env_finalise


contains
    subroutine parse_command_line
        integer            :: i
        character(len=512) :: arg

        niter = 10
        i = 0
        offset = 0
        nfiles = 0
        comm_type = 'p2p'
        l_subcomm = .false.
        ncfname = ''

        do
            call get_command_argument(i, arg)
            if (len_trim(arg) == 0) then
                exit
            endif

            if (arg == '--basename') then
                i = i + 1
                call get_command_argument(i, arg)
                basename = trim(arg)
            else if (arg == '--niter') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') niter
            else if (arg == '--offset') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') offset
            else if (arg == '--nfiles') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') nfiles
            else if (arg == '--size_factor') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') parcel%size_factor
            else if (arg == '--comm-type') then
                i = i + 1
                call get_command_argument(i, arg)
                comm_type = trim(arg)
#ifdef ENABLE_COARRAY
                call mpi_print("WARNING: Ignoring 'comm_type' argument. Coarray is enabled.")
#endif
            else if (arg == '--subcomm') then
                l_subcomm = .true.
            else if (arg == '--ncfname') then
                i = i + 1
                call get_command_argument(i, arg)
                ncfname = trim(arg)
            else if (arg == '--help') then
                if (world%rank == world%root) then
                    print *, "./benchmark_read ",                               &
                             "--basename [basename] ",                          &
                             "--niter [int] ",                                  &
                             "--offset [int] ",                                 &
                             "--nfiles [int] ",                                 &
                             "--subcomm (optional, disabled for shmem) ",       &
                             "--comm-type [p2p, rma, shmem] ",                  &
                             "--ncfname [string]",                              &
                             "--size_factor [float]"
                endif
                call mpi_stop
            endif
            i = i+1
        enddo

        if (ncfname == '') then
            call mpi_stop("No netCDF output file name provided.")
        endif

#ifdef ENABLE_COARRAY
        comm_type = 'caf'
#endif

        if (world%rank == world%root) then
            print *, "basename", basename
            print *, "offset", offset
            print *, "number of files", nfiles
            print *, "niter", niter
            print *, "size_factor", parcel%size_factor
            print *, "enabled subcommunicator", l_subcomm
            print *, "comm type: " // comm_type
            print *, "netCDF output file name: " // trim(ncfname)
        endif

    end subroutine parse_command_line

end program benchmark_read
