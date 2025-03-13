! =============================================================================
!               Timer module implemented according to PMPIC
!                     https://github.com/EPCCed/pmpic
! =============================================================================
module mpi_timer
    use mpi_environment
    use mpi_collectives
    use mpi_f08, only : MPI_Wtime
    implicit none

    type timer_type
        character(len=32)   :: name
        integer             :: handle = -1
        double precision    :: wall_time
        logical             :: running
        integer             :: n_calls
        double precision    :: start_time, end_time
        double precision    :: mean_time
        double precision    :: min_time
        double precision    :: max_time
    end type timer_type

    type(timer_type), allocatable :: timings(:)

    integer :: n_timers = 0

    private :: n_timers, get_statistics

contains

    subroutine register_timer(name, handle)
        character(*), intent(in)    :: name
        integer,      intent(inout) :: handle

        if (handle /= -1) then
            ! Timer already registered.
            return
        endif

        n_timers = n_timers + 1

        handle = n_timers

        if (n_timers > size(timings)) then
            call mpi_stop("Timer number exceeds allocated size.")
        endif


        timings(handle)%name = name
        timings(handle)%handle = handle
        timings(handle)%wall_time = 0.0d0
        timings(handle)%running = .false.
        timings(handle)%n_calls = 0
    end subroutine register_timer

    subroutine start_timer(handle)
        integer, intent(in) :: handle

        if (handle == -1) then
            call mpi_exit_on_error("Trying to start unregistered timer.")
        endif

        if (timings(handle)%running) then
            return
        endif

        timings(handle)%n_calls = timings(handle)%n_calls + 1

        timings(handle)%running = .true.

        timings(handle)%start_time = MPI_Wtime()

    end subroutine start_timer

    subroutine stop_timer(handle)
        integer, intent(in) :: handle

        if (handle == -1) then
            call mpi_exit_on_error("Trying to stop unregistered timer.")
        endif

        if (.not. timings(handle)%running) then
            return
        endif

        timings(handle)%running = .false.

        timings(handle)%end_time = MPI_Wtime()

        timings(handle)%wall_time = timings(handle)%wall_time &
                                  + timings(handle)%end_time  &
                                  - timings(handle)%start_time
    end subroutine stop_timer

    subroutine write_timings(fname)
        character(*), intent(in)     :: fname
        character(len=len(fname)+12) :: csv_timer_file
        character(len=len(fname)+11) :: csv_ncall_file
        logical                      :: l_exist = .false.
        character(len=3)             :: status = 'new'
        integer                      :: i

        timings(1:n_timers)%max_time = get_statistics(MPI_MAX)

        ! we need to take the maximum number because of the subcommunicator in the
        ! parcel nearest algorithm
        call mpi_blocking_reduce(timings(1:n_timers)%n_calls, MPI_MAX, world)

        if (.not. world%rank == world%root) then
            return
        endif

        csv_timer_file = trim(fname) // '-timings.csv'

        inquire(file=csv_timer_file, exist=l_exist)

        if (l_exist) then
            status = 'old'
        endif

        open(unit=1234, file=csv_timer_file, status=status, position='append')

        if (.not. l_exist) then
            do i = 1, n_timers-1
                write(1234, '(a)', advance='no') trim(timings(i)%name) // ','
            enddo
            write(1234, '(a)', advance='yes') trim(timings(n_timers)%name)
        endif

        do i = 1, n_timers-1
            write(1234, '(es24.16e3,a1)', advance='no') timings(i)%max_time, ","
        enddo
        write(1234, '(es24.16e3)') timings(n_timers)%max_time

        close(1234)


        csv_ncall_file = trim(fname) // '-ncalls.csv'

        inquire(file=csv_ncall_file, exist=l_exist)

        if (l_exist) then
            status = 'old'
        endif

        open(unit=1235, file=csv_ncall_file, status=status, position='append')

        if (.not. l_exist) then
            do i = 1, n_timers-1
                write(1235, '(a)', advance='no') trim(timings(i)%name) // ','
            enddo
            write(1235, '(a)', advance='yes') trim(timings(n_timers)%name)
        endif

        do i = 1, n_timers-1
            write(1235, "(i0,a1)", advance='no') timings(i)%n_calls, ","
        enddo
        write(1235, "(i0)") timings(n_timers)%n_calls

        close(1235)


    end subroutine write_timings


    function get_statistics(op) result(buffer)
        type(MPI_Op), intent(in) :: op
        double precision         :: buffer(1:n_timers)

        buffer = timings(1:n_timers)%wall_time

        if (world%rank == world%root) then
            call MPI_Reduce(MPI_IN_PLACE,           &
                            buffer(1:n_timers),     &
                            n_timers,               &
                            MPI_DOUBLE_PRECISION,   &
                            op,                     &
                            world%root,             &
                            world%comm,             &
                            world%err)

        else
            call MPI_Reduce(buffer(1:n_timers),     &
                            buffer(1:n_timers),     &
                            n_timers,               &
                            MPI_DOUBLE_PRECISION,   &
                            op,                     &
                            world%root,             &
                            world%comm,             &
                            world%err)
        endif

    end function get_statistics

end module mpi_timer
