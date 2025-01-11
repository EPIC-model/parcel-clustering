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
    use utils, only : total_timer              &
                    , register_timer           &
                    , register_all_timers      &
                    , print_timer              &
                    , start_timer              &
                    , stop_timer
    use parcel_merging
    use parcel_netcdf
    use netcdf_utils
    use netcdf_reader
    use iomanip, only : zfill
    use parcel_nearest_p2p_graph, only : p2p_graph_t
    use parcel_nearest_rma_graph, only : rma_graph_t
    use parcel_nearest_shmem_graph, only : shmem_graph_t
    implicit none

    integer          :: allreduce_timer = -1
    character(64)    :: basename
    character(512)   :: fname
    integer          :: ncid, n, m, niter
    integer          :: ncells(3), offset, nfiles
    character(len=9) :: graph_type ! OpenSHMEM, MPI RMA, MPI P2P
    logical          :: l_subcomm

    call mpi_env_initialise

    call register_all_timers
    call register_timer('MPI allreduce', allreduce_timer)

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

    select case(graph_type)
        case ('MPI P2P')
            allocate(p2p_graph_t :: tree)
        case ('MPI RMA')
            allocate(rma_graph_t :: tree)
        case ('OpenSHMEM')
            allocate(shmem_graph_t :: tree)
        case default
            allocate(p2p_graph_t :: tree)
    end select

    call tree%initialise(max_num_parcels, l_subcomm)

    call tree%register_timer

    call start_timer(total_timer)

    do n = 0, niter - 1

        m = mod(n, nfiles) + offset

        fname = trim(basename) // '_' // zfill(m) // '_parcels.nc'
        if (world%rank == world%root) then
           print *, "Read:", trim(fname)
        endif
        call read_netcdf_parcels(fname)

        parcels%total_num = 0

        call start_timer(allreduce_timer)
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

        call stop_timer(allreduce_timer)

        call parcel_merge

        parcels%total_num = 0
        call start_timer(allreduce_timer)
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

        call stop_timer(allreduce_timer)
    enddo

    call stop_timer(total_timer)

    call parcels%deallocate

    call print_timer

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
        graph_type = 'MPI P2P'
        l_subcomm = .false.

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
            else if (arg == '--graph-type') then
                i = i + 1
                call get_command_argument(i, arg)
                graph_type = trim(arg)
            else if (arg == '--subcomm') then
                l_subcomm = .true.
            else if (arg == '--help') then
                if (world%rank == world%root) then
                    print *, "./benchmark_read ",                               &
                             "--basename [basename] ",                          &
                             "--niter [int] ",                                  &
                             "--offset [int] ",                                 &
                             "--nfiles [int] ",                                 &
                             "--subcomm (optional, disabled for OpenhSHMEM) ",  &
                             "--graph-type [MPI P2P, MPI RMA, OpenSHMEM] ",     &
                             "--size_facctor [float]"
                endif
                call mpi_stop
            endif
            i = i+1
        enddo

        if (world%rank == world%root) then
            print *, "basename", basename
            print *, "offset", offset
            print *, "number of files", nfiles
            print *, "niter", niter
            print *, "size_factor", parcel%size_factor
            print *, "enabled subcommunicator", l_subcomm
            print *, "graph type: " // graph_type
        endif

    end subroutine parse_command_line

end program benchmark_read
