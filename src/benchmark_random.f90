program benchmark_random
    use parcel_container
    use options, only : parcel
    use parameters, only : update_parameters, lower, extent, nx, ny, nz, max_num_parcels, vmin
    use parcel_init, only : parcel_default
    use parcel_merging
    use parcel_nearest, only : tree
    use parcel_mpi, only : parcel_communicate
    use mpi_environment
    use mpi_layout
    use mpi_timer, only : write_timings
    use mpi_datatypes, only : MPI_INTEGER_64BIT
    use mpi_ops, only : MPI_SUM_64BIT
    use mpi_utils, only : mpi_stop
#ifdef ENABLE_COARRAY
    use mpi_utils, only : mpi_print
#endif
    use utils, only : register_all_timers      &
                    , setup_parcels            &
                    , init_rng
#ifndef ENABLE_COARRAY
    use parcel_nearest_p2p_graph, only : p2p_graph_t
    use parcel_nearest_rma_graph, only : rma_graph_t
#ifdef ENABLE_SHMEM
    use parcel_nearest_shmem_graph, only : shmem_graph_t
#endif
#endif
    implicit none

    integer             :: k, niter, seed
    integer(kind=int64) :: n_small_parcels, n_remaining_parcels
    double precision    :: lx, ly, lz, small_parcel_fraction
    logical             :: l_shuffle, l_variable_nppc, l_subcomm
    character(len=512)  :: csvfname
    character(len=5)    :: comm_type ! shmem, rma, p2p or caf
    character(len=1)    :: snum
    integer(kind=int64) :: buf(9) ! size(n_way_parcel_mergers) = 7; +1 (n_parcel_merges); +1 (n_big_close)

    call mpi_env_initialise


    call parse_command_line

    lower  = (/zero, zero, zero/)
    extent = (/lx, ly, lz/)

    parcel%lambda_max = 4.0d0

    call mpi_layout_init(lower, extent, nx, ny, nz)

    call update_parameters

    call parcels%allocate(max_num_parcels)

    call init_rng(seed)

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

    call register_all_timers

    call tree%register_timer

    do k = 1, niter

        ! -------------------------------------------------------------
        ! Set up the parcel configuration:
        call setup_parcels(small_parcel_fraction, l_shuffle, l_variable_nppc)

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

        call parcel_communicate(parcels)

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
!     n_big_close = buf(2)
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
        integer            :: i
        character(len=512) :: arg

        nx = 32
        ny = 32
        nz = 32
        lx = 128.0d0
        ly = 128.0d0
        lz = 128.0d0
        niter = 1
        small_parcel_fraction = 0.5d0
        parcel%n_per_cell = 40
        parcel%min_vratio = 40.0d0
        parcel%size_factor = 1.25d0
        l_shuffle = .false.
        l_variable_nppc = .false.
        l_subcomm = .false.
        comm_type = 'p2p'
        seed = 42
        csvfname = ''


        i = 1
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
            else if (arg == '--niter') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') niter
            else if (arg == '--nppc') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') parcel%n_per_cell
            else if (arg == '--min-vratio') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') parcel%min_vratio
            else if (arg == '--size-factor') then
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
            else if (arg == '--small-parcel-fraction') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(f16.0)') small_parcel_fraction
            else if (arg == '--shuffle') then
                l_shuffle = .true.
            else if (arg == '--variable-nppc') then
                l_variable_nppc = .true.
            else if (arg == '--subcomm') then
                l_subcomm = .true.
            else if (arg == '--comm-type') then
                i = i + 1
                call get_command_argument(i, arg)
                comm_type = trim(arg)
#ifdef ENABLE_COARRAY
                call mpi_print("WARNING: Ignoring 'comm_type' argument. Coarray is enabled.")
#endif
            else if (arg == '--seed') then
                i = i + 1
                call get_command_argument(i, arg)
                read(arg,'(i6)') seed
            else if (arg == '--csvfname') then
                i = i + 1
                call get_command_argument(i, arg)
                csvfname = trim(arg)
            else if (arg == '--help') then
                if (world%rank == world%root) then
                    print *, "./benchmark_random ",                              &
                             "--nx [int] --ny [int] --nz [int] ",                &
                             "--lx [float] --ly [float] --lz [float] ",          &
                             "--niter [int] --nppc [int] ",                      &
                             "--min-vratio [float] --size-factor [float] ",      &
                             "--small-parcel-fration [float] ",                  &
                             "--shuffle (optional) --variable-nppc (optional) ", &
                             "--subcomm (optional, disabled for shmem) ",        &
                             "--seed [int] ",                                    &
                             "--csvfname [string]",                              &
                             "--comm-type [p2p, rma, shmem]"
                endif
                call mpi_stop
            else
                call mpi_stop("Unknown input argument.")
            endif
            i = i+1
        end do

        if (csvfname == '') then
            call mpi_stop("No timing output file name provided.")
        endif

#ifdef ENABLE_COARRAY
        comm_type = "caf"
#endif

        if (world%rank == world%root) then
            print *, "nx", nx
            print *, "ny", ny
            print *, "nz", nz
            print *, "lx", lx
            print *, "ly", ly
            print *, "lz", lz
            print *, "niter", niter
            print *, "nppc", parcel%n_per_cell
            print *, "min_vratio", parcel%min_vratio
            print *, "size_factor", parcel%size_factor
            print *, "fraction of small parcels (%)", small_parcel_fraction * 100.d0
            print *, "shuffle parcels", l_shuffle
            print *, "seed", seed
            print *, "enabled subcommunicator", l_subcomm
            print *, "variable number of parcels/cell:", l_variable_nppc
            print *, "comm type: " // comm_type
            print *, "ASCII file name: " // trim(csvfname)
        endif

    end subroutine parse_command_line

end program benchmark_random
