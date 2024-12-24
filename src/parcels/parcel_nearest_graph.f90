module parcel_nearest_graph
    use mpi_layout
    use mpi_utils
    use parcel_mpi, only : get_parcel_id_buffer_ptr     &
                         , deallocate_parcel_id_buffers
    use datatypes, only : intlog_pair_t
    use mpi_datatypes, only : MPI_INTEGER_LOGICAL_ARRAY
    implicit none

    private

    type, abstract :: graph_t

    contains
        procedure(graph_initialise),     deferred :: initialise
        procedure(graph_finalise),       deferred :: finalise
        procedure(graph_reset),          deferred :: reset
        procedure(graph_resolve),        deferred :: resolve
        procedure(graph_register_timer), deferred :: register_timer
    end type

    interface
        subroutine graph_initialise(this, num)
            import :: graph_t
            class(graph_t), intent(inout) :: this
            integer,        intent(in)    :: num
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
        subroutine graph_resolve(this, mpi_comm, isma, iclo, rclo, n_local_small)
            use mpi_environment, only : communicator
            import :: graph_t
            class(graph_t),     intent(inout) :: this
            type(communicator), intent(inout) :: mpi_comm
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

end module parcel_nearest_graph
