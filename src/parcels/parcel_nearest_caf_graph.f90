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

        logical :: l_caf_allocated = .false.

        ! Mapping of neighbouring ranks between MPI Cartesian topology and CAF
        integer :: cart2caf(8)
        integer :: me ! this PE

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

        procedure, private :: init_cart2caf
        procedure, private :: get_pe

    end type

    public :: caf_graph_t

contains

    subroutine caf_graph_initialise(this, num, l_subcomm)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: num
        logical,            intent(in)    :: l_subcomm

        if (this%l_caf_allocated) then
            return
        endif

        if (.not. l_mpi_layout_initialised) then
            call mpi_stop("Error: The Cartesian communicator not yet initialised.")
        endif

        this%l_caf_allocated = .true.

        ! Ensure we use all MPI ranks because we need to call
        ! sync
        if (l_subcomm) then
            call mpi_print("Ignoring the request to use a subcommunicator.")
            call mpi_print("We only support a global communicator with Coarray Fortran.")
        endif
        this%l_enabled_subcomm = .false.

        !--------------------------------------------------
        ! Check mapping beteen MPI and CAF

        call this%init_cart2caf

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

        this%l_merged(:)[this%me] = .false.
        this%l_leaf(:)[this%me] = .false.
        this%l_available(:)[this%me] = .false.

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
                if (.not. this%l_merged(is)[this%me]) then
                    ic = iclo(m)
                    rc = rclo(m)
                    this%l_leaf(is)[this%me] = .true.
                    call this%put_avail(rc, ic, .true.)
                endif
            enddo

            ! This barrier is necessary!
            call this%barrier

            ! determine leaf parcels
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)[this%me]) then
                    ic = iclo(m)
                    rc = rclo(m)
                    call this%put_leaf(rc, ic, .false.)
                endif
            enddo

            ! We must synchronise all MPI processes here to ensure all MPI processes
            ! have done theirRMA operations as we modify the windows again.
            call this%barrier

            ! filter out parcels that are "unavailable" for merging
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)[this%me]) then
                    if (.not. this%l_leaf(is)[this%me]) then
                        ic = iclo(m)
                        rc = rclo(m)
                        call this%put_avail(rc, ic, .false.)
                    endif
                endif
            enddo

            ! This sync is necessary as CAF processes access their l_available
            ! array which may be modified in the loop above. In order to make sure all
            ! CAF processes have finished above loop, we need this barrier.
            call this%barrier

            ! identify mergers in this iteration
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)[this%me]) then
                    ic = iclo(m)
                    rc = rclo(m)

                    l_helper = this%get_avail(rc, ic)

                    if (this%l_leaf(is)[this%me] .and. l_helper) then
                        l_continue_iteration = .true. ! merger means continue iteration
                        this%l_merged(is)[this%me] = .true.

                        call this%put_merged(rc, ic, .true.)
                    endif
                endif
            enddo

            call start_timer(this%allreduce_timer)
            ! Performance improvement: We actually only need to synchronize with neighbours
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               l_continue_iteration,    &
                               1,                       &
                               MPI_LOGICAL,             &
                               MPI_LOR,                 &
                               this%comm%comm,          &
                               this%comm%err)
            call stop_timer(this%allreduce_timer)
            call mpi_check_for_error(this%comm, &
                "in MPI_Allreduce of caf_graph_t::resolve_tree.")
        enddo

        ! No barrier necessary because of the blocking MPI_Allreduce that acts like
        ! a barrier!

        ! Second stage, related to dual links
        do m = 1, n_local_small
            is = isma(m)

            if (.not. this%l_merged(is)[this%me]) then
                if (this%l_leaf(is)[this%me]) then ! set in last iteration of stage 1
                    ic = iclo(m)
                    rc = rclo(m)
                    call this%put_avail(rc, ic, .true.)
                endif
            endif
        enddo

        ! This barrier is necessary as we modifiy l_available above and need it below
        call this%barrier

        ! Second stage
        do m = 1, n_local_small
            is = isma(m)
            ic = iclo(m)
            rc = rclo(m)
            l_do_merge(m) = .false.
            l_isolated_dual_link(m) = .false.

            if (this%l_merged(is)[this%me] .and. this%l_leaf(is)[this%me]) then
                ! previously identified mergers: keep
                l_do_merge(m) = .true.
                !----------------------------------------------------------
                ! begin of sanity check
                ! After first stage mergers parcel cannot be both initiator
                ! and receiver in stage 1
                l_helper = this%get_leaf(rc, ic)

                if (l_helper) then
                    call mpi_exit_on_error(&
                        'in caf_graph_t::resolve_tree: First stage error')
                endif

                ! end of sanity check
                !----------------------------------------------------------

            elseif (.not. this%l_merged(is)[this%me]) then
                if (this%l_leaf(is)[this%me]) then
                    ! links from leafs
                    l_do_merge(m) = .true.
                elseif (.not. this%l_available(is)[this%me]) then
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
                            this%l_available(is)[this%me] = .true.
                        endif
                    endif
                endif
            endif
        enddo


        ! This barrier is necessary.
        call this%barrier

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

        call register_timer('resolve graphs', this%resolve_timer)
        call register_timer('MPI allreduce', this%allreduce_timer)
        call register_timer('Coarray put', this%put_timer)
        call register_timer('Coarray get', this%get_timer)
        call register_timer('Coarray sync', this%sync_timer)

    end subroutine caf_graph_register_timer

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_avail(this, rank, ic, val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val
        integer                           :: pe

        if (rank == cart%rank) then
            this%l_available(ic)[this%me] = val
        else
            call start_timer(this%put_timer)
            pe = this%get_pe(rank)
            this%l_available(ic)[pe] = val
            call stop_timer(this%put_timer)
        endif

    end subroutine put_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_leaf(this, rank, ic, val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val
        integer                           :: pe

        if (rank == cart%rank) then
            this%l_leaf(ic)[this%me] = val
        else
            call start_timer(this%put_timer)
            pe = this%get_pe(rank)
            this%l_leaf(ic)[pe] = val
            call stop_timer(this%put_timer)
        endif

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val
        integer                           :: pe

        if (rank == cart%rank) then
            this%l_merged(ic)[this%me] = val
        else
            call start_timer(this%put_timer)
            pe = this%get_pe(rank)
            this%l_merged(ic)[pe] = val
            call stop_timer(this%put_timer)
        endif

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_avail(this, rank, ic) result(val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val
        integer                           :: pe

        if (rank == cart%rank) then
            val = this%l_available(ic)[this%me]
        else
            call start_timer(this%get_timer)
            pe = this%get_pe(rank)
            val = this%l_available(ic)[pe]
            call stop_timer(this%get_timer)
        endif

    end function get_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val
        integer                           :: pe

        if (rank == cart%rank) then
            val = this%l_leaf(ic)[this%me]
        else
            call start_timer(this%get_timer)
            pe = this%get_pe(rank)
            val = this%l_leaf(ic)[pe]
            call stop_timer(this%get_timer)
        endif

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(this, rank, ic) result(val)
        class(caf_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val
        integer                           :: pe

        if (rank == cart%rank) then
            val = this%l_merged(ic)[this%me]
        else
            call start_timer(this%get_timer)
            pe = this%get_pe(rank)
            val = this%l_merged(ic)[pe]
            call stop_timer(this%get_timer)
        endif

    end function get_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine barrier(this)
        class(caf_graph_t), intent(inout) :: this

        call start_timer(this%sync_timer)
        sync images([this%cart2caf])
        call stop_timer(this%sync_timer)

    end subroutine barrier

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_cart2caf(this)
        class(caf_graph_t), intent(inout) :: this
        type(MPI_Request)                 :: requests(8)
        type(MPI_Status)                  :: statuses(8)
        integer                           :: n

        ! *this* CAF PE (= processing element)
        this%me = this_image()

        ! check if MPI comm world rank == CAF rank
        ! Note: MPI rank starts at 0; CAF rank starts at 1
        if (world%rank + 1 /= this%me) then
            call mpi_exit_on_error("MPI rank and CAF do not agree.")
        endif

        this%cart2caf = -1

        do n = 1, 8
            call MPI_Isend(this%me,                 &
                           1,                       &
                           MPI_INTEGER,             &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of caf_graph_t::init_cart2caf.")
        enddo

        do n = 1, 8
            call MPI_Recv(this%cart2caf(n),       &
                          1,                      &
                          MPI_INTEGER,            &
                          neighbours(n)%rank,     &
                          RECV_NEIGHBOUR_TAG(n),  &
                          cart%comm,              &
                          statuses(n),            &
                          cart%err)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of caf_graph_t::init_cart2caf.")

        if (any(this%cart2caf == -1)) then
            call mpi_exit_on_error("in caf_graph_t::init_cart2caf: Not all MPI ranks finished.")
        endif

    end subroutine init_cart2caf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_pe(this, rank) result(pe)
        class(caf_graph_t), intent(in) :: this
        integer,            intent(in) :: rank
        integer                        :: n
        integer                        :: pe

        n = get_neighbour_from_rank(rank)
        pe = this%cart2caf(n)

    end function get_pe

end module parcel_nearest_caf_graph
