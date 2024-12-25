module parcel_nearest_rma_graph
    use mpi_layout
    use mpi_utils
    use mpi_timer, only : start_timer       &
                        , stop_timer        &
                        , register_timer
    use parcel_nearest_graph, only : graph_t
    use iso_c_binding, only : c_ptr, c_f_pointer
    implicit none

    private

    type, extends(graph_t) :: rma_graph_t

        ! Logicals used to determine which mergers are executed
        ! Integers above could be reused for this, but this would
        ! make the algorithm less readable
        logical, pointer :: l_leaf(:)
        logical, pointer :: l_available(:)
        logical, pointer :: l_merged(:)    ! indicates parcels merged in first stage

        integer :: resolve_timer
        integer :: allreduce_timer
        integer :: rma_put_timer
        integer :: rma_get_timer
        integer :: sync_timer

        type(MPI_Win) :: win_merged, win_avail, win_leaf
        logical       :: l_win_allocated

    contains

        procedure :: initialise => rma_graph_initialise
        procedure :: finalise   => rma_graph_finalise
        procedure :: reset      => rma_graph_reset
        procedure :: resolve    => rma_graph_resolve
        procedure :: register_timer => rma_graph_register_timer

        procedure, private :: put_avail
        procedure, private :: put_leaf
        procedure, private :: put_merged

        procedure, private :: get_avail
        procedure, private :: get_leaf
        procedure, private :: get_merged

        procedure, private :: barrier

    end type

    public :: rma_graph_t

contains

    subroutine rma_graph_initialise(this, num)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: num
        integer (KIND=MPI_ADDRESS_KIND)   :: win_size
        logical                           :: l_bytes
        integer                           :: disp_unit
        integer (KIND=MPI_ADDRESS_KIND)   :: long_type
        type(c_ptr)                       :: buf_ptr

        if (this%l_win_allocated) then
            return
        endif

        if (.not. l_mpi_layout_initialised) then
            call mpi_stop("Error: The Cartesian communicator not yet initialised.")
        endif

        this%l_win_allocated = .true.

        call MPI_Sizeof(l_bytes, disp_unit, cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Sizeof of rma_graph_t::initialise.")

        ! size of RMA window in bytes
        long_type = disp_unit
        win_size = long_type * num

        if (win_size < 0) then
            call mpi_stop("Error: Integer overflow. Unable to allocate MPI RMA windows.")
        endif

        ! allocate window win_leaf and memory for l_leaf
        call MPI_Win_allocate(win_size,         &
                              disp_unit,        &
                              MPI_INFO_NULL,    &
                              cart%comm,        &
                              buf_ptr,          &
                              this%win_leaf,    &
                              cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Win_allocate of rma_graph_t::initialise.")

        call c_f_pointer(buf_ptr, this%l_leaf, [num])


        ! allocate window win_avail and memory for l_available
        call MPI_Win_allocate(win_size,         &
                              disp_unit,        &
                              MPI_INFO_NULL,    &
                              cart%comm,        &
                              buf_ptr,          &
                              this%win_avail,   &
                              cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Win_allocate of rma_graph_t::initialise.")

        call c_f_pointer(buf_ptr, this%l_available, [num])

        ! allocate window win_merged and memory for l_merged
        call MPI_Win_allocate(win_size,         &
                              disp_unit,        &
                              MPI_INFO_NULL,    &
                              cart%comm,        &
                              buf_ptr,          &
                              this%win_merged,  &
                              cart%err)

        call mpi_check_for_error(cart, &
            "in MPI_Win_allocate of rma_graph_t::initialise.")

        call c_f_pointer(buf_ptr, this%l_merged, [num])

        call mpi_check_rma_window_model(this%win_avail)
        call mpi_check_rma_window_model(this%win_merged)
        call mpi_check_rma_window_model(this%win_leaf)

    end subroutine rma_graph_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine rma_graph_finalise(this)
        class(rma_graph_t), intent(inout) :: this

        if (.not. this%l_win_allocated) then
            return
        endif

        call MPI_Win_free(this%win_leaf, cart%err)
        call mpi_check_for_error(cart, &
                "in MPI_Win_free of rma_graph_t::finalise.")

        call MPI_Win_free(this%win_avail, cart%err)
        call mpi_check_for_error(cart, &
                "in MPI_Win_free of rma_graph_t::finalise.")

        call MPI_Win_free(this%win_merged, cart%err)
        call mpi_check_for_error(cart, &
                "in MPI_Win_free of rma_graph_t::finalise.")

    end subroutine rma_graph_finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine rma_graph_reset(this)
        class(rma_graph_t), intent(inout) :: this

        this%l_merged = .false.
        this%l_leaf = .false.
        this%l_available = .false.

    end subroutine rma_graph_reset

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! https://github.com/mpi-forum/mpi-forum-historic/issues/413
    ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node294.htm
    ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node279.htm
    subroutine rma_graph_resolve(this, mpi_comm, isma, iclo, rclo, n_local_small)
        class(rma_graph_t), intent(inout) :: this
        type(communicator), intent(inout) :: mpi_comm
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
!             call MPI_Barrier(mpi_comm%comm, mpi_comm%err)
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
!             call MPI_Barrier(mpi_comm%comm, mpi_comm%err)
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

            ! This MPI_Barrier is necessary as MPI processes access their l_available
            ! array which may be modified in the loop above. In order to make sure all
            ! MPI ranks have finished above loop, we need this barrier.
            call start_timer(this%sync_timer)
            call this%barrier
!             call MPI_Barrier(mpi_comm%comm, mpi_comm%err)
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
            ! Performance improvement: We actually only need to synchronize with neighbours
            call MPI_Allreduce(MPI_IN_PLACE,            &
                               l_continue_iteration,    &
                               1,                       &
                               MPI_LOGICAL,             &
                               MPI_LOR,                 &
                               mpi_comm%comm,           &
                               mpi_comm%err)
            call stop_timer(this%allreduce_timer)
            call mpi_check_for_error(mpi_comm, &
                "in MPI_Allreduce of parcel_nearest::resolve_tree.")
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
!         call MPI_Barrier(mpi_comm%comm, mpi_comm%err)
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
                        'in parcel_nearest::resolve_tree: First stage error')
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
!         call MPI_Barrier(mpi_comm%comm, mpi_comm%err)
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
    end subroutine rma_graph_resolve

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine rma_graph_register_timer(this)
        class(rma_graph_t), intent(inout) :: this

        call register_timer('graph resolve', this%resolve_timer)
        call register_timer('MPI graph allreduce', this%allreduce_timer)
        call register_timer('MPI RMA put', this%rma_put_timer)
        call register_timer('MPI RMA get', this%rma_get_timer)
        call register_timer('MPI graph sync', this%sync_timer)

    end subroutine rma_graph_register_timer

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_avail(this, rank, ic, val)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val
        integer(KIND=MPI_ADDRESS_KIND)    :: offset

        if (rank == cart%rank) then
            this%l_available(ic) = val
        else
            call start_timer(this%rma_put_timer)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this%win_avail, cart%err)
            !     MPI_Put(origin_addr, origin_count, origin_datatype, target_rank,
            !         target_disp, target_count, target_datatype, win, ierror)
            !     TYPE(*), DIMENSION(..), INTENT(IN), ASYNCHRONOUS :: origin_addr
            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
            !     TYPE(MPI_Win), INTENT(IN) :: win
            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
            offset = ic - 1 ! starts at 0
            call MPI_Put(val,              &
                         1,                &
                         MPI_LOGICAL,      &
                         rank,             &
                         offset,           &
                         1,                &
                         MPI_LOGICAL,      &
                         this%win_avail,   &
                         cart%err)
            call mpi_check_for_error(cart, &
                "in MPI_Put of parcel_nearest::resolve_tree.")

            ! After MPI_Win_unlock, the RMA operation is completed at the origin and target.
            call MPI_Win_unlock(rank, this%win_avail, cart%err)
            call stop_timer(this%rma_put_timer)
        endif

    end subroutine put_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_leaf(this, rank, ic, val)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val
        integer(KIND=MPI_ADDRESS_KIND)    :: offset

        if (rank == cart%rank) then
            this%l_leaf(ic) = .false.
        else
            call start_timer(this%rma_put_timer)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this%win_leaf, cart%err)
            offset = ic - 1
            call MPI_Put(val,              &
                         1,                &
                         MPI_LOGICAL,      &
                         rank,             &
                         offset,           &
                         1,                &
                         MPI_LOGICAL,      &
                         this%win_leaf,    &
                         cart%err)
            call mpi_check_for_error(cart, &
                "in MPI_Put of parcel_nearest::resolve_tree.")
            call MPI_Win_unlock(rank, this%win_leaf, cart%err)
            call stop_timer(this%rma_put_timer)
        endif

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical,            intent(in)    :: val
        integer(KIND=MPI_ADDRESS_KIND)    :: offset

        if (rank == cart%rank) then
            this%l_merged(ic) = .true.
        else
            call start_timer(this%rma_put_timer)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this%win_merged, cart%err)

            offset = ic - 1
            call MPI_Put(val,              &
                         1,                &
                         MPI_LOGICAL,      &
                         rank,             &
                         offset,           &
                         1,                &
                         MPI_LOGICAL,      &
                         this%win_merged,  &
                         cart%err)
            call mpi_check_for_error(cart, &
                "in MPI_Put of parcel_nearest::resolve_tree.")
            call MPI_Win_unlock(rank, this%win_merged, cart%err)
            call stop_timer(this%rma_put_timer)
        endif

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_avail(this, rank, ic) result(val)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val
        integer(KIND=MPI_ADDRESS_KIND)    :: offset

        if (rank == cart%rank) then
            val = this%l_available(ic)
        else
            call start_timer(this%rma_get_timer)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this%win_avail, cart%err)
            ! Note: The OpenMPI specification says that processes must be on the same node
            !       in order MPI_Get to work. However, I tested a simple MPI_Get on Archer2
            !       between nodes with 1 rank per node. It works! It may therefore only be
            !       a limitation of OpenMPI.
            !     MPI_Get(origin_addr, origin_count, origin_datatype, target_rank,
            !         target_disp, target_count, target_datatype, win, ierror)
            !     TYPE(*), DIMENSION(..), ASYNCHRONOUS :: origin_addr
            !     INTEGER, INTENT(IN) :: origin_count, target_rank, target_count
            !     TYPE(MPI_Datatype), INTENT(IN) :: origin_datatype, target_datatype
            !     INTEGER(KIND=MPI_ADDRESS_KIND), INTENT(IN) :: target_disp
            !     TYPE(MPI_Win), INTENT(IN) :: win
            !     INTEGER, OPTIONAL, INTENT(OUT) :: ierror
            offset = ic - 1
            call MPI_Get(val,              &
                         1,                &
                         MPI_LOGICAL,      &
                         rank,             &
                         offset,           &
                         1,                &
                         MPI_LOGICAL,      &
                         this%win_avail,   &
                          cart%err)
            call mpi_check_for_error(cart, &
                    "in MPI_Get of parcel_nearest::resolve_tree.")

            call MPI_Win_unlock(rank, this%win_avail, cart%err)
            call stop_timer(this%rma_get_timer)
        endif

    end function get_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val
        integer(KIND=MPI_ADDRESS_KIND)    :: offset

        if (rank == cart%rank) then
            val = this%l_leaf(ic)
        else
            call start_timer(this%rma_get_timer)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this%win_leaf, cart%err)
            offset = ic - 1
            call MPI_Get(val,              &
                         1,                &
                         MPI_LOGICAL,      &
                         rank,             &
                         offset,           &
                         1,                &
                         MPI_LOGICAL,      &
                         this%win_leaf,    &
                         cart%err)
            call mpi_check_for_error(cart, &
                "in MPI_Get of parcel_nearest::resolve_tree.")

            call MPI_Win_unlock(rank, this%win_leaf, cart%err)
            call stop_timer(this%rma_get_timer)
        endif

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(this, rank, ic) result(val)
        class(rma_graph_t), intent(inout) :: this
        integer,            intent(in)    :: rank
        integer,            intent(in)    :: ic
        logical                           :: val
        integer(KIND=MPI_ADDRESS_KIND)    :: offset

        if (rank == cart%rank) then
            val = this%l_merged(ic)
        else
            call start_timer(this%rma_get_timer)
            call MPI_Win_lock(MPI_LOCK_SHARED, rank, 0, this%win_merged, cart%err)
            offset = ic - 1
            call MPI_Get(val,              &
                         1,                &
                         MPI_LOGICAL,      &
                         rank,             &
                         offset,           &
                         1,                &
                         MPI_LOGICAL,      &
                         this%win_merged,  &
                         cart%err)
            call mpi_check_for_error(cart, &
                "in MPI_Get of parcel_nearest::resolve_tree.")

            call MPI_Win_unlock(rank, this%win_merged, cart%err)
            call stop_timer(this%rma_get_timer)
        endif

    end function get_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine barrier(this)
        class(rma_graph_t), intent(inout) :: this
        type(MPI_Request)                 :: requests(8)
        type(MPI_Status)                  :: statuses(8)
        integer                           :: n
        logical                           :: l_send, l_recv(8)

        l_send = .true.

        !----------------------------------------------------------------------
        ! Send from remote to owning rank and sync data at owning rank
        do n = 1, 8
            call MPI_Isend(l_send,                  &
                           1,                       &
                           MPI_LOGICAL,             &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of rma_graph_t::barrier.")
        enddo

        do n = 1, 8
            call MPI_Recv(l_recv(n),                &
                          1,                        &
                          MPI_LOGICAL,              &
                          neighbours(n)%rank,       &
                          RECV_NEIGHBOUR_TAG(n),    &
                          cart%comm,                &
                          statuses(n),              &
                          cart%err)
        enddo

        call MPI_Waitall(8,                 &
                        requests,           &
                        statuses,           &
                        cart%err)

        call mpi_check_for_error(cart, &
                                "in MPI_Waitall of rma_graph_t::barrier.")

        if (.not. all(l_recv)) then
            call mpi_exit_on_error("in rma_graph_t::barrier: Not all MPI ranks finished.")
        endif

    end subroutine barrier

end module parcel_nearest_rma_graph
