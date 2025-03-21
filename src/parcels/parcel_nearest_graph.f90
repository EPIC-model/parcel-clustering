module parcel_nearest_graph
    use mpi_environment
    use mpi_layout, only : cart, neighbours
    use mpi_utils, only : mpi_check_for_error   &
                        , mpi_exit_on_error
    use mpi_timer, only : start_timer       &
                        , stop_timer        &
                        , register_timer
    implicit none

    private

    type, abstract :: graph_t

        logical :: l_enabled_subcomm = .false.

        type(communicator) :: comm

#ifdef ENABLE_COARRAY
        ! Because the extending type caf_graph_t has a coarray component, the
        ! parent type ‘graph_t’ must have one as well
        ! (see also, https://fortran-lang.discourse.group/t/extensible-derived-type-that-contains-a-coarray/1691
        ! [16 Jan 2025]
        logical, codimension[:], allocatable :: l_dummy
#endif

        integer :: resolve_timer = -1
        integer :: allreduce_timer = -1
        integer :: put_timer = -1
        integer :: get_timer = -1
        integer :: sync_timer = -1
        integer :: comm_timer = -1

        ! each rank can have at most 8 neighbours
        logical :: l_valid_comm_neighbour(8)

    contains
        procedure(graph_initialise),     deferred :: initialise
        procedure(graph_finalise),       deferred :: finalise
        procedure(graph_reset),          deferred :: reset
        procedure(graph_resolve),        deferred :: resolve
        procedure(graph_register_timer), deferred :: register_timer

        procedure :: create_comm
        procedure :: register_common_timers
    end type

    interface
        subroutine graph_initialise(this, num, l_subcomm)
            import :: graph_t
            class(graph_t), intent(inout) :: this
            integer,        intent(in)    :: num
            logical,        intent(in)    :: l_subcomm
        end subroutine graph_initialise
    end interface

    interface
        subroutine graph_finalise(this)
            import :: graph_t
            class(graph_t), intent(inout) :: this
        end subroutine graph_finalise
    end interface

    interface
        subroutine graph_reset(this)
            import :: graph_t
            class(graph_t), intent(inout) :: this
        end subroutine graph_reset
    end interface

    interface
        subroutine graph_resolve(this, isma, iclo, rclo, n_local_small)
            use mpi_environment, only : communicator
            import :: graph_t
            class(graph_t),     intent(inout) :: this
            integer,            intent(inout) :: isma(0:)
            integer,            intent(inout) :: iclo(:)
            integer,            intent(inout) :: rclo(:)
            integer,            intent(inout) :: n_local_small
        end subroutine graph_resolve
    end interface

    interface
        subroutine graph_register_timer(this)
            import :: graph_t
            class(graph_t), intent(inout) :: this
        end subroutine graph_register_timer
    end interface

    public :: graph_t

contains

    subroutine create_comm(this, l_exclude)
        class(graph_t), intent(inout) :: this
        logical,        intent(in)    :: l_exclude
        integer                       :: color, n
        type(MPI_Request)             :: requests(8)
        type(MPI_Status)              :: statuses(8)
        logical                       :: l_valid

        call start_timer(this%comm_timer)

        if (.not. this%l_enabled_subcomm) then
            call MPI_Comm_dup(cart%comm, this%comm%comm, this%comm%err)
            call mpi_check_for_error(this%comm, &
                    "in MPI_Comm_dup of graph_t::create_comm.")

            call MPI_Comm_size(this%comm%comm, this%comm%size, this%comm%err)
            call mpi_check_for_error(this%comm, &
                    "in MPI_Comm_size of graph_t::create_comm.")

            if (this%comm%size /= cart%size) then
                call mpi_exit_on_error("MPI Cartesian size not identical.")
            endif

            call MPI_Comm_rank(this%comm%comm, this%comm%rank, this%comm%err)
            call mpi_check_for_error(this%comm, &
                    "in MPI_Comm_rank of graph_t::create_comm.")

            if (this%comm%rank /= cart%rank) then
                call mpi_exit_on_error("MPI Cartesian rank not identical.")
            endif

            this%comm%root = cart%root

            this%l_valid_comm_neighbour = .true.

            return
        endif

        ! Ensure the communicator is freed first.
        if (this%comm%comm /= MPI_COMM_NULL) then
            call MPI_Comm_free(this%comm%comm, this%comm%err)
            call mpi_check_for_error(this%comm, &
                    "in MPI_Comm_free of graph_t::create_comm.")
        endif

        ! Each MPI process must know if it is part of the this%communicator or not.
        ! All MPI ranks that have small parcels or received small parcels from neighbouring
        ! MPI ranks must be part of the communicator.
        color = MPI_UNDEFINED
        if (.not. l_exclude) then
            color = 0  ! any non-negative number is fine
        endif

        call MPI_Comm_split(comm=cart%comm,          &
                            color=color,             &
                            key=cart%rank,           &  ! key controls the ordering of the processes
                            newcomm=this%comm%comm,  &
                            ierror=cart%err)

        if (this%comm%comm /= MPI_COMM_NULL) then
            ! The following two calls are not necessary, but we do for good practice.
            call MPI_Comm_size(this%comm%comm, this%comm%size, this%comm%err)
            call mpi_check_for_error(this%comm, &
                    "in MPI_Comm_size of graph_t::create_comm.")
            call MPI_Comm_rank(this%comm%comm, this%comm%rank, this%comm%err)
            call mpi_check_for_error(this%comm, &
                    "in MPI_Comm_rank of graph_t::create_comm.")
            this%comm%root = 0
        endif

        ! Figure out which neighbour is part of sub-communicator:
        l_valid = (this%comm%comm /= MPI_COMM_NULL)

        do n = 1, 8
            call MPI_Isend(l_valid,                 &
                           1,                       &
                           MPI_LOGICAL,             &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of graph_t::barrier.")
        enddo

        do n = 1, 8
            call MPI_Recv(this%l_valid_comm_neighbour(n), &
                          1,                              &
                          MPI_LOGICAL,                    &
                          neighbours(n)%rank,             &
                          RECV_NEIGHBOUR_TAG(n),          &
                          cart%comm,                      &
                          statuses(n),                    &
                          cart%err)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Waitall of graph_t::create_comm.")

        call stop_timer(this%comm_timer)

    end subroutine create_comm

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine register_common_timers(this, label)
        class(graph_t), intent(inout) :: this
        character(*),   intent(in)    :: label

        call register_timer('create comm', this%comm_timer)
        call register_timer('resolve graphs', this%resolve_timer)
        call register_timer('MPI allreduce', this%allreduce_timer)
        call register_timer(label // ' put', this%put_timer)
        call register_timer(label // ' get', this%get_timer)
        call register_timer(label // ' sync', this%sync_timer)

    end subroutine register_common_timers

end module parcel_nearest_graph
