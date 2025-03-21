program benchmark_read
    use constants, only : pi, zero, one, f12, f23, twopi
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels, vmin
    use parcel_init, only : parcel_default
    use parcel_nearest, only : tree
    use parcel_mpi, only : parcel_communicate
    use mpi_environment
    use mpi_timer, only : write_timings
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
    implicit none

    character(427)      :: dirname
    character(64)       :: basename
    character(512)      :: csvfname, ncfname
    integer             :: ncid, n, m, niter, k
    integer             :: ncells(3), offset, nfiles
    integer(kind=int64) :: n_small_parcels, n_remaining_parcels
    character(len=5)    :: comm_type ! p2p, rma, shmem or caf
    logical             :: l_subcomm
    character(len=1)    :: snum
    integer(kind=int64) :: buf(9)

    call mpi_env_initialise

    call register_all_timers

    parcel%lambda_max = 4.0d0
    parcel%min_vratio = 20.0d0
    parcel%size_factor = 1.0d0

    call parse_command_line

    ncfname = trim(dirname) // trim(basename) // '_' // zfill(offset) // '_parcels.nc'
    call open_netcdf_file(trim(ncfname), NF90_NOWRITE, ncid)

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

        ncfname = trim(dirname) // trim(basename) // '_' // zfill(m) // '_parcels.nc'
        if (world%rank == world%root) then
           print *, "Read: ", trim(ncfname)
        endif
        call read_netcdf_parcels(ncfname)

        n_small_parcels = count(parcels%volume(1:parcels%local_num) < vmin, kind=int64)

        call mpi_blocking_reduce(n_small_parcels, MPI_SUM, world)

        parcels%total_num = 0

        call MPI_Allreduce(parcels%local_num, &
                           parcels%total_num, &
                           1,                 &
                           MPI_INTEGER_64BIT, &
                           MPI_SUM_64BIT,     &
                           world%comm,        &
                           world%err)

        if (world%rank == world%root) then
            n_remaining_parcels = parcels%total_num
            print *, "Number of parcels before merging: ", parcels%total_num
            print '(a,f8.4,a)', " Fraction of small parcels:                    ", &
                                n_small_parcels / dble(parcels%total_num) * 100.0d0, "%"
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
            print *, "Number of parcels after merging:  ", parcels%total_num
            print '(a,f8.4,a)', " Fraction of merged parcels:                   ", &
                (n_remaining_parcels - parcels%total_num) / dble(n_remaining_parcels) * 100.d0, "%"
        endif
    enddo

    call parcels%deallocate

    buf(1) = n_parcel_merges
    buf(2) = n_big_close
    buf(3:9) = n_way_parcel_mergers

    call mpi_blocking_reduce(buf, MPI_SUM, world)

    n_parcel_merges = buf(1)

    n_way_parcel_mergers = buf(3:9)

    if (world%rank == world%root) then
        print *, "Number of MPI ranks:        ", world%size
        print *, "Total number of merges:     ", n_parcel_merges
        print *, "Number of close big parcels:", buf(2) !n_big_close
        do k = 1, 7
            write(snum, fmt='(I1)')  k+1
            print *, "Number of " // snum // "-way mergers:    ", n_way_parcel_mergers(k)
        enddo
    endif

    call write_timings(trim(csvfname))

    call tree%finalise

    call mpi_env_finalise


contains
    subroutine parse_command_line
        integer            :: i, dirlen
        character(len=512) :: arg

        niter = 10
        offset = 0
        nfiles = 0
        comm_type = 'p2p'
        l_subcomm = .false.
        csvfname = ''

        i = 1
        do
            call get_command_argument(i, arg)
            if (len_trim(arg) == 0) then
                exit
            endif

            if (arg == '--dirname') then
                i = i + 1
                call get_command_argument(i, arg)
                dirname = trim(arg)
            else if (arg == '--ncbasename') then
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
            else if (arg == '--size-factor') then
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
            else if (arg == '--csvfname') then
                i = i + 1
                call get_command_argument(i, arg)
                csvfname = trim(arg)
            else if (arg == '--help') then
                if (world%rank == world%root) then
                    print *, "./benchmark_read ",                               &
                             "--dirname [directory] ",                          &
                             "--ncbasename [basename] ",                        &
                             "--niter [int] ",                                  &
                             "--offset [int] ",                                 &
                             "--nfiles [int] ",                                 &
                             "--subcomm (optional, disabled for shmem) ",       &
                             "--comm-type [p2p, rma, shmem] ",                  &
                             "--csvfname [string]",                             &
                             "--size-factor [float]"
                endif
                call mpi_stop
            else
                call mpi_stop("Unknown input argument.")
            endif
            i = i+1
        enddo

        if (csvfname == '') then
            call mpi_stop("No timing output file name provided.")
        endif

#ifdef ENABLE_COARRAY
        comm_type = 'caf'
#endif

        dirlen = len(trim(dirname))

        if (dirname(dirlen:dirlen) /= '/') then
            dirname = trim(dirname) // '/'
        endif


        if (world%rank == world%root) then
            print *, "dirname                 ", trim(dirname)
            print *, "basename                ", trim(basename)
            print *, "offset                  ", offset
            print *, "number of files         ", nfiles
            print *, "niter                   ", niter
            print *, "size-factor             ", parcel%size_factor
            print *, "enabled subcommunicator ", l_subcomm
            print *, "comm-type:              " // comm_type
            print *, "ASCII file name:        " // trim(csvfname)
        endif

    end subroutine parse_command_line

end program benchmark_read
