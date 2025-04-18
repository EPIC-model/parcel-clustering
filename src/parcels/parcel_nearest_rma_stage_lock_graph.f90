module parcel_nearest_rma_stage_lock_graph
    use mpi_layout
    use mpi_utils
    use mpi_timer, only : start_timer       &
                        , stop_timer        &
                        , register_timer
    use parcel_nearest_rma_graph, only : rma_graph_t
    use iso_c_binding, only : c_ptr, c_f_pointer
    implicit none

    private

    type, extends(rma_graph_t) :: rma_stage_lock_graph_t

    contains

        procedure :: resolve    => rma_stage_lock_graph_resolve
        procedure :: register_timer => rma_stage_lock_graph_register_timer

        procedure :: put_avail
        procedure :: put_leaf
        procedure :: put_merged

        procedure :: get_avail
        procedure :: get_leaf
        procedure :: get_merged

    end type

    public :: rma_stage_lock_graph_t

contains

    ! https://github.com/mpi-forum/mpi-forum-historic/issues/413
    ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node294.htm
    ! https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node279.htm
    subroutine rma_stage_lock_graph_resolve(this, isma, iclo, rclo, n_local_small)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(inout) :: isma(0:)
        integer,                       intent(inout) :: iclo(:)
        integer,                       intent(inout) :: rclo(:)
        integer,                       intent(inout) :: n_local_small
        integer                                      :: ic, rc, is, m, j
        logical                                      :: l_helper
        logical                                      :: l_continue_iteration, l_do_merge(n_local_small)
        logical                                      :: l_isolated_dual_link(n_local_small)

        call start_timer(this%resolve_timer)

        ! First, iterative, stage
        l_continue_iteration = .true.

        do while (l_continue_iteration)
            l_continue_iteration = .false.
            ! reset relevant properties for candidate mergers

            call MPI_Win_lock_all(0, this%win_avail)

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

            call MPI_Win_unlock_all(this%win_avail)

            ! This barrier is necessary!
            call this%barrier

            call MPI_Win_lock_all(0, this%win_leaf)

            ! determine leaf parcels
            do m = 1, n_local_small
                is = isma(m)

                if (.not. this%l_merged(is)) then
                    ic = iclo(m)
                    rc = rclo(m)
                    call this%put_leaf(rc, ic, .false.)
                endif
            enddo

            call MPI_Win_unlock_all(this%win_leaf)

            ! We must synchronise all MPI processes here to ensure all MPI processes
            ! have done theirRMA operations as we modify the windows again.
            call this%barrier

            call MPI_Win_lock_all(0, this%win_avail)

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

            call MPI_Win_unlock_all(this%win_avail)

            ! This MPI_Barrier is necessary as MPI processes access their l_available
            ! array which may be modified in the loop above. In order to make sure all
            ! MPI ranks have finished above loop, we need this barrier.
            call this%barrier

            call MPI_Win_lock_all(0, this%win_avail)
            call MPI_Win_lock_all(0, this%win_merged)

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

            call MPI_Win_unlock_all(this%win_merged)
            call MPI_Win_unlock_all(this%win_avail)

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
                "in MPI_Allreduce of parcel_nearest::resolve_tree.")
        enddo

        ! No barrier necessary because of the blocking MPI_Allreduce that acts like
        ! a barrier!

        call MPI_Win_lock_all(0, this%win_avail)

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

        call MPI_Win_unlock_all(this%win_avail)

        ! This barrier is necessary as we modifiy l_available above and need it below.
        call this%barrier

        call MPI_Win_lock_all(0, this%win_leaf)
        call MPI_Win_lock_all(0, this%win_avail)

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

        call MPI_Win_unlock_all(this%win_avail)
        call MPI_Win_unlock_all(this%win_leaf)

        ! This barrier is necessary.
        call this%barrier

        call MPI_Win_lock_all(0, this%win_avail)

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

        call MPI_Win_unlock_all(this%win_avail)

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
    end subroutine rma_stage_lock_graph_resolve

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine rma_stage_lock_graph_register_timer(this)
        class(rma_stage_lock_graph_t), intent(inout) :: this

        call this%register_common_timers(label='MPI RMA (stage lock)')

    end subroutine rma_stage_lock_graph_register_timer

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_avail(this, rank, ic, val)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(in)    :: rank
        integer,                       intent(in)    :: ic
        logical,                       intent(in)    :: val
        integer(KIND=MPI_ADDRESS_KIND)               :: offset

        if (rank == cart%rank) then
            this%l_available(ic) = val
        else
            call start_timer(this%put_timer)
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

            ! Complete RMA operation at the origin and the target
            call MPI_Win_flush(rank, this%win_avail, cart%err)

            call stop_timer(this%put_timer)
        endif

    end subroutine put_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_leaf(this, rank, ic, val)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(in)    :: rank
        integer,                       intent(in)    :: ic
        logical,                       intent(in)    :: val
        integer(KIND=MPI_ADDRESS_KIND)               :: offset

        if (rank == cart%rank) then
            this%l_leaf(ic) = val
        else
            call start_timer(this%put_timer)
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

            call MPI_Win_flush(rank, this%win_leaf, cart%err)

            call stop_timer(this%put_timer)
        endif

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(in)    :: rank
        integer,                       intent(in)    :: ic
        logical,                       intent(in)    :: val
        integer(KIND=MPI_ADDRESS_KIND)               :: offset

        if (rank == cart%rank) then
            this%l_merged(ic) = val
        else
            call start_timer(this%put_timer)

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

            call MPI_Win_flush(rank, this%win_merged, cart%err)

            call stop_timer(this%put_timer)
        endif

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_avail(this, rank, ic) result(val)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(in)    :: rank
        integer,                       intent(in)    :: ic
        logical                                      :: val
        integer(KIND=MPI_ADDRESS_KIND)               :: offset


        if (rank == cart%rank) then
            val = this%l_available(ic)
        else
            call start_timer(this%get_timer)
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

            call MPI_Win_flush(rank, this%win_avail, cart%err)

            call stop_timer(this%get_timer)
        endif

    end function get_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(in)    :: rank
        integer,                       intent(in)    :: ic
        logical                                      :: val
        integer(KIND=MPI_ADDRESS_KIND)               :: offset

        if (rank == cart%rank) then
            val = this%l_leaf(ic)
        else
            call start_timer(this%get_timer)
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

            call MPI_Win_flush(rank, this%win_leaf, cart%err)

            call stop_timer(this%get_timer)
        endif

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(this, rank, ic) result(val)
        class(rma_stage_lock_graph_t), intent(inout) :: this
        integer,                       intent(in)    :: rank
        integer,                       intent(in)    :: ic
        logical                                      :: val
        integer(KIND=MPI_ADDRESS_KIND)               :: offset

        if (rank == cart%rank) then
            val = this%l_merged(ic)
        else
            call start_timer(this%get_timer)
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

            call MPI_Win_flush(rank, this%win_merged, cart%err)

            call stop_timer(this%get_timer)
        endif

    end function get_merged

end module parcel_nearest_rma_stage_lock_graph
