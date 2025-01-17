module parcel_nearest_graph
    use mpi_environment
    use mpi_layout, only : cart
    use mpi_utils, only : mpi_check_for_error   &
                        , mpi_exit_on_error
    implicit none

    private

    type, abstract :: graph_t

        logical :: l_enabled_subcomm = .false.

        type(communicator) :: comm

        integer :: resolve_timer = -1
        integer :: allreduce_timer = -1
        integer :: put_timer = -1
        integer :: get_timer = -1
        integer :: sync_timer = -1

    contains
        procedure(graph_initialise),     deferred :: initialise
        procedure(graph_finalise),       deferred :: finalise
        procedure(graph_reset),          deferred :: reset
        procedure(graph_resolve),        deferred :: resolve
        procedure(graph_register_timer), deferred :: register_timer

        procedure :: create_comm
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

    subroutine create_comm(this, l_include)
        class(graph_t), intent(inout) :: this
        logical,        intent(in)    :: l_include
        integer                       :: color

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
        if (.not. l_include) then
            color = 0  ! any non-negative number is fine
        endif

        call MPI_Comm_split(comm=cart%comm,         &
                            color=color,            &
                            key=cart%rank,          &  ! key controls the ordering of the processes
                            newcomm=this%comm%comm,   &
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

    end subroutine create_comm

end module parcel_nearest_graph
