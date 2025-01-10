program benchmark_verify
    use constants, only : pi, zero, one, f12, f23, twopi
    use parcels_mod, only : parcels
    use options, only : parcel, output
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels
    use parcel_init, only : parcel_default
    use parcel_mpi, only : parcel_communicate
    use mpi_environment
    use mpi_layout
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_ops, only : MPI_SUM_64BIT
    use mpi_utils, only : mpi_stop
    use utils, only : epic_timer               &
                    , register_timer           &
                    , register_all_timers      &
                    , print_timer              &
                    , start_timer              &
                    , stop_timer               &
                    , setup_parcels            &
                    , init_rng
    use parcel_netcdf
    use parcel_nearest_p2p_graph, only : p2p_graph_t
    use parcel_nearest_rma_graph, only : rma_graph_t
    use parcel_nearest_shmem_graph, only : shmem_graph_t
    use parcel_nearest, only : tree
    implicit none

    integer              :: allreduce_timer, seed
    double precision     :: lx, ly, lz
    logical              :: l_setup, l_shuffle, l_subcomm
    character(len=9)     :: graph_type ! OpenSHMEM, MPI RMA, MPI P2P

    call mpi_env_initialise

    call register_all_timers
    call register_timer('MPI allreduce', allreduce_timer)

    call parse_command_line

    lower  = (/zero, zero, zero/)
    extent = (/lx, ly, lz/)

    parcel%lambda_max = 4.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call parcels%allocate(max_num_parcels)

    call init_rng(seed)

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

    call start_timer(epic_timer)

    ! -------------------------------------------------------------
    ! Set up the parcel configuration:
    if (l_setup) then
        call setup_parcels(xlen=lx, ylen=ly, zlen=lz, l_shuffle=l_shuffle, l_variable_nppc=.false.)

        call start_timer(allreduce_timer)
        parcels%total_num = 0
        call MPI_Allreduce(parcels%local_num, &
                           parcels%total_num, &
                           1,                 &
                           MPI_INTEGER_64BIT, &
                           MPI_SUM_64BIT,     &
                           world%comm,        &
                           world%err)
        call stop_timer(allreduce_timer)

        output%parcel_list(1)  = 'volume'
        output%parcel_list(2)  = 'x_position'
        output%parcel_list(3)  = 'y_position'
        output%parcel_list(4)  = 'z_position'
        output%parcel_list(5)  = 'buoyancy'
        output%parcel_list(6)  = 'x_vorticity'
        output%parcel_list(7)  = 'y_vorticity'
        output%parcel_list(8)  = 'z_vorticity'
        output%parcel_list(9)  = 'B11'
        output%parcel_list(10) = 'B12'
        output%parcel_list(11) = 'B13'
        output%parcel_list(12) = 'B22'
        output%parcel_list(13) = 'B23'
        call create_netcdf_parcel_file('initial', .true., .false.)
        call write_netcdf_parcels(t = 0.0d0)
    else
        call read_netcdf_parcels('initial_0000000001_parcels.nc')
    endif


    parcels%total_num = 0

    call start_timer(allreduce_timer)
    call MPI_Allreduce(parcels%local_num, &
                       parcels%total_num, &
                       1,                 &
                       MPI_INTEGER_64BIT, &
                       MPI_SUM_64BIT,     &
                       world%comm,        &
                       world%err)

    call stop_timer(allreduce_timer)

    if (world%size == 1) then
        call serial_merge
    else
        call parallel_merge
    endif

    parcels%total_num = 0
    call start_timer(allreduce_timer)
    call MPI_Allreduce(parcels%local_num, &
                       parcels%total_num, &
                       1,                 &
                       MPI_INTEGER_64BIT, &
                       MPI_SUM_64BIT,     &
                       world%comm,        &
                       world%err)
    call stop_timer(allreduce_timer)


    if (world%size == 1) then
        call create_netcdf_parcel_file('serial_final', .true., .false.)
    else
        call create_netcdf_parcel_file('parallel_final', .true., .false.)
    endif

    call write_netcdf_parcels(t = 0.0d0)

    call stop_timer(epic_timer)

    call parcels%deallocate

    call print_timer

    call tree%finalise

    call mpi_env_finalise

contains
    subroutine parse_command_line
        integer            :: i
        character(len=512) :: arg

        nx = 32
        ny = 32
        nz = 32
        lx = 1.0d0
        ly = 1.0d0
        lz = 1.0d0
        parcel%n_per_cell = 40
        parcel%min_vratio = 40.0d0
        parcel%size_factor = 1.25d0
        l_setup = .false.
        l_shuffle = .false.
        l_subcomm = .false.
        graph_type = 'MPI P2P'
        seed = 42


        i = 0
        do
            call get_command_argument(i, arg)
            if (len_trim(arg) == 0) then
                exit
            endif

            if (arg == '--nx') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') nx
            else if (arg == '--ny') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') ny
            else if (arg == '--nz') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') nz
            else if (arg == '--n_per_cell') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') parcel%n_per_cell
            else if (arg == '--min_vratio') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') parcel%min_vratio
            else if (arg == '--size_factor') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') parcel%size_factor
            else if (arg == '--lx') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') lx
            else if (arg == '--ly') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') ly
            else if (arg == '--lz') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') lz
            else if (arg == '--setup-parcels') then
                l_setup = .true.
            else if (arg == '--shuffle') then
                l_shuffle = .true.
            else if (arg == '--subcomm') then
                l_subcomm = .true.
            else if (arg == '--graph-type') then
                i = i + 1
                call get_command_argument(i, arg)
                graph_type = trim(arg)
            else if (arg == '--seed') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') seed
            else if (arg == '--help') then
                if (world%rank == world%root) then
                    print *, "./benchmark_verify --nx [int] --ny [int] --nz [int] ",  &
                             "--n_per_cell [int] ",                                   &
                             "--min_vratio [float] ",                                 &
                             "--shuffle (optional) ",                                 &
                             "--seed [int] ",                                         &
                             "--subcomm (optional, disabled for OpenhSHMEM) ",        &
                             "--graph-type [MPI P2P, MPI RMA, OpenSHMEM]"
                endif
                call mpi_stop
            endif
            i = i+1
        end do

        if (world%rank == world%root) then
            print *, "nx", nx
            print *, "ny", ny
            print *, "nz", nz
            print *, "lx", lx
            print *, "ly", ly
            print *, "lz", lz
            print *, "n_per_cell", parcel%n_per_cell
            print *, "min_vratio", parcel%min_vratio
            print *, "size_factor", parcel%size_factor
            print *, "shuffle parcels", l_shuffle
            print *, "seed", seed
            print *, "enabled subcommunicator", l_subcomm
            print *, "graph type: " // graph_type
        endif

    end subroutine parse_command_line

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine serial_merge
        use parcel_merge_serial
        use parcel_nearest_serial, only : merge_nearest_timer       &
                                        , merge_tree_resolve_timer

        call register_timer('parcel merge', merge_timer)
        call register_timer('merge nearest', merge_nearest_timer)
        call register_timer('merge tree resolve', merge_tree_resolve_timer)

        call merge_parcels

    end subroutine serial_merge

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine parallel_merge
        use parcel_merging

        call register_timer('parcel merge', merge_timer)
        call register_timer('merge nearest', merge_nearest_timer)

        call tree%reset

        call parcel_merge

    end subroutine parallel_merge

end program benchmark_verify
