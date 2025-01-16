module parcel_nearest_caf_graph
    use mpi_layout
    use mpi_utils
    use mpi_timer, only : start_timer       &
                        , stop_timer        &
                        , register_timer
    use parcel_nearest_graph, only : graph_t
    use mpi_environment, only : l_ignore_mpi_finalize
    implicit none

    private

    type, extends(graph_t) :: caf_graph_t

        private

        ! Logicals used to determine which mergers are executed
        ! Integers above could be reused for this, but this would
        ! make the algorithm less readable
        logical, dimension(:), codimension[:], allocatable :: l_leaf
        logical, dimension(:), codimension[:], allocatable :: l_available
        logical, dimension(:), codimension[:], allocatable :: l_merged    ! indicates parcels merged in first stage

        integer :: resolve_timer = -1
        integer :: allreduce_timer = -1
        integer :: caf_put_timer = -1
        integer :: caf_get_timer = -1
        integer :: sync_timer = -1

        logical :: l_caf_allocated = .false.

    contains

        procedure :: initialise => caf_graph_initialise
        procedure :: finalise   => caf_graph_finalise
        procedure :: reset      => caf_graph_reset
        procedure :: resolve    => caf_graph_resolve
        procedure :: register_timer => caf_graph_register_timer

        procedure, private :: put_avail
        procedure, private :: put_leaf
        procedure, private :: put_merged

        procedure, private :: get_avail
        procedure, private :: get_leaf
        procedure, private :: get_merged

        procedure, private :: barrier

    end type

    public :: caf_graph_t

contains

    subroutine caf_graph_initialise(this, num, l_subcomm)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: num
        logical,            intent(in)    :: l_subcomm
!         integer                           :: error

        if (this%l_caf_allocated) then
            return
        endif

        if (.not. l_mpi_layout_initialised) then
            call mpi_stop("Error: The Cartesian communicator not yet initialised.")
        endif

        this%l_caf_allocated = .true.

        ! Ensure we use all MPI ranks because we need to call
        ! shmem_barrier_all
        if (l_subcomm) then
            call mpi_print("Ignoring the request to use a subcommunicator.")
            call mpi_print("We only support a global communicator with OpenSHMEM.")
        endif
        this%l_enabled_subcomm = .false.

        !--------------------------------------------------

        allocate (this%l_available(num)[*])
        allocate (this%l_leaf(num)[*])
        allocate (this%l_merged(num)[*])

        call this%reset

    end subroutine caf_graph_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine caf_graph_finalise(this)
        class(caf_graph_t), intent(inout) :: this
        type(c_ptr)                       :: buf_ptr

        if (.not. this%l_caf_allocated) then
            return
        endif

        this%l_caf_allocated = .false.

        sync images(*)

        l_ignore_mpi_finalize = .true.

    end subroutine caf_graph_finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine caf_graph_reset(this)
        class(caf_graph_t), intent(inout) :: this

        this%l_merged = .false.
        this%l_leaf = .false.
        this%l_available = .false.

    end subroutine caf_graph_reset

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine caf_graph_resolve(this, isma, iclo, rclo, n_local_small)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(inout) :: isma(0:)
        integer,            intent(inout) :: iclo(:)
        integer,            intent(inout) :: rclo(:)
        integer,            intent(inout) :: n_local_small
        integer                           :: ic, rc, is, m, j
        logical                           :: l_helper
        logical                           :: l_continue_iteration, l_do_merge(n_local_small)
        logical                           :: l_isolated_dual_link(n_local_small)

        call start_timer(this%resolve_timer)

        ! First, iterative, stage
        l_continue_iteration = .true.

        do while (l_continue_iteration)
            l_continue_iteration = .false.
            ! reset relevant properties for candidate mergers

            do m = 1, n_local_small
                is = isma(m)
                ! only consider links that still may be merging
                ! reset relevant properties
                if (.not. this%l_merged(is)) then
                    ic = iclo(m)
                    rc = rclo(m)
                    this%l_leaf(is) = .true.
                    call this%put_avail(rc, ic, .true.)
                endif
            enddo

            ! This barrier is necessary!
            call start_timer(this%sync_timer)
            call this%barrier
            call stop_timer(this%sync_timer)

            ! determine leaf parcels
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)) then
                    ic = iclo(m)
                    rc = rclo(m)
                    call this%put_leaf(rc, ic, .false.)
                endif
            enddo

            ! We must synchronise all MPI processes here to ensure all MPI processes
            ! have done theirRMA operations as we modify the windows again.
            call start_timer(this%sync_timer)
            call this%barrier
            call stop_timer(this%sync_timer)

            ! filter out parcels that are "unavailable" for merging
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)) then
                    if (.not. this%l_leaf(is)) then
                        ic = iclo(m)
                        rc = rclo(m)
                        call this%put_avail(rc, ic, .false.)
                    endif
                endif
            enddo

            ! This sync is necessary as SHMEM processes access their l_available
            ! array which may be modified in the loop above. In order to make sure all
            ! SHMEM processes have finished above loop, we need this barrier.
            call start_timer(this%sync_timer)
            call this%barrier
            call stop_timer(this%sync_timer)


            ! identify mergers in this iteration
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)) then
                    ic = iclo(m)
                    rc = rclo(m)

                    l_helper = this%get_avail(rc, ic)

                    if (this%l_leaf(is) .and. l_helper) then
                        l_continue_iteration = .true. ! merger means continue iteration
                        this%l_merged(is) = .true.

                        call this%put_merged(rc, ic, .true.)
                    endif
                endif
            enddo

            call start_timer(this%allreduce_timer)
            ! Perfoshmemnce improvement: We actually only need to synchronize with neighbours
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               l_continue_iteration,    &
                               1,                       &
                               MPI_LOGICAL,             &
                               MPI_LOR,                 &
                               this%comm%comm,          &
                               this%comm%err)
            call stop_timer(this%allreduce_timer)
            call mpi_check_for_error(this%comm, &
                "in MPI_Allreduce of shmem_graph_t::resolve_tree.")
        enddo

        ! No barrier necessary because of the blocking MPI_Allreduce that acts like
        ! a barrier!

        ! Second stage, related to dual links
        do m = 1, n_local_small
            is = isma(m)

            if (.not. this%l_merged(is)) then
                if (this%l_leaf(is)) then ! set in last iteration of stage 1
                    ic = iclo(m)
                    rc = rclo(m)
                    call this%put_avail(rc, ic, .true.)
                endif
            endif
        enddo

        ! This barrier is necessary as we modifiy l_available above and need it below.
        call start_timer(this%sync_timer)
        call this%barrier
        call stop_timer(this%sync_timer)

        ! Second stage
        do m = 1, n_local_small
            is = isma(m)
            ic = iclo(m)
            rc = rclo(m)
            l_do_merge(m) = .false.
            l_isolated_dual_link(m) = .false.

            if (this%l_merged(is) .and. this%l_leaf(is)) then
                ! previously identified mergers: keep
                l_do_merge(m) = .true.
                !----------------------------------------------------------
                ! begin of sanity check
                ! After first stage mergers parcel cannot be both initiator
                ! and receiver in stage 1
                l_helper = this%get_leaf(rc, ic)

                if (l_helper) then
                    call mpi_exit_on_error(&
                        'in shmem_graph_t::resolve_tree: First stage error')
                endif

                ! end of sanity check
                !----------------------------------------------------------

            elseif (.not. this%l_merged(is)) then
                if (this%l_leaf(is)) then
                    ! links from leafs
                    l_do_merge(m) = .true.
                elseif (.not. this%l_available(is)) then
                    ! Above means parcels that have been made 'available' do not keep outgoing links

                    l_helper = this%get_avail(rc, ic)

                    if (l_helper) then
                        ! merge this parcel into ic along with the leaf parcels
                        l_do_merge(m) = .true.

                    else
                        l_isolated_dual_link(m) = .true.
                        ! isolated dual link
                        ! Don't keep current link
                        ! But make small parcel available so other parcel can merge with it
                        ! THIS NEEDS THINKING ABOUT A PARALLEL IMPLEMENTATION
                        ! This could be based on the other parcel being outside the domain
                        ! And a "processor order"
                        if (cart%rank <= rc) then
                            ! The MPI rank with lower number makes its parcel
                            ! available.
                            this%l_available(is) = .true.
                        endif
                    endif
                endif
            endif
        enddo


        ! This barrier is necessary.
        call start_timer(this%sync_timer)
        call this%barrier
        call stop_timer(this%sync_timer)

        !------------------------------------------------------
        do m = 1, n_local_small
            is = isma(m)
            ic = iclo(m)
            rc = rclo(m)

            if ((l_do_merge(m) .eqv. .false.) .and. l_isolated_dual_link(m)) then
                ! isolated dual link

                l_helper = this%get_avail(rc, ic)

                if (l_helper) then
                    ! merge this parcel into ic along with the leaf parcels
                    l_do_merge(m) = .true.
                !else
                !   ! Dual link is resolved on other rank
                endif
            endif
            !------------------------------------------------------
        enddo

        j = 0
        do m = 1, n_local_small
            is = isma(m)
            ic = iclo(m)
            rc = rclo(m)
            if (l_do_merge(m)) then
                j = j + 1
                isma(j) = is
                iclo(j) = ic
                rclo(j) = rc
            endif
        enddo
        n_local_small = j

        call stop_timer(this%resolve_timer)
    end subroutine caf_graph_resolve

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine caf_graph_register_timer(this)
        class(caf_graph_t), intent(inout) :: this

        call register_timer('graph resolve', this%resolve_timer)
        call register_timer('MPI graph allreduce', this%allreduce_timer)
        call register_timer('Coarray put', this%caf_put_timer)
        call register_timer('Coarray get', this%caf_get_timer)
        call register_timer('Coarray sync', this%sync_timer)

    end subroutine caf_graph_register_timer

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_avail(this, rank, ic, val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val

        call start_timer(this%caf_put_timer)
        this%l_available(ic) = val
        call stop_timer(this%caf_put_timer)

    end subroutine put_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_leaf(this, rank, ic, val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val

        call start_timer(this%caf_put_timer)
        this%l_leaf(ic) = val
        call stop_timer(this%caf_put_timer)

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val

        call start_timer(this%caf_put_timer)
        this%l_merged(ic) = val
        call stop_timer(this%caf_put_timer)

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_avail(this, rank, ic) result(val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val

        call start_timer(this%caf_get_timer)
        val = this%l_available(ic)
        call stop_timer(this%caf_get_timer)

    end function get_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val

        call start_timer(this%caf_get_timer)
        val = this%l_leaf(ic)
        call stop_timer(this%caf_get_timer)

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(this, rank, ic) result(val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val

        call start_timer(this%caf_get_timer)
        val = this%l_merged(ic)
        call stop_timer(this%caf_get_timer)

    end function get_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine barrier(this)
        class(caf_graph_t), intent(inout) :: this

        sync images(*)

    end subroutine barrier

end module parcel_nearest_caf_graph
