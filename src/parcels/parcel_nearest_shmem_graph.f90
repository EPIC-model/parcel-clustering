module shmem
    implicit none

    interface
        subroutine shmem_init() bind(C, name="shmem_init")
            use, intrinsic :: iso_c_binding
        end subroutine shmem_init

        subroutine shmem_finalize() bind(C, name="shmem_finalize")
            use, intrinsic :: iso_c_binding
        end subroutine shmem_finalize

        type(c_ptr) function c_shmem_malloc(n) bind(C, name="shmem_malloc")
            use, intrinsic :: iso_c_binding
            integer(kind=c_size_t), value :: n
        end function c_shmem_malloc

        subroutine c_shmem_free(ptr) bind(C, name="shmem_free")
            use, intrinsic :: iso_c_binding
            type(c_ptr), value :: ptr
        end subroutine c_shmem_free

        integer(kind=c_int) function shmem_pe_accessible(pe) bind(C, name="shmem_pe_accessible")
            use, intrinsic :: iso_c_binding
            integer(kind=c_int), value :: pe
        end function shmem_pe_accessible

        integer(kind=c_int) function shmem_addr_accessible(addr, pe) bind(C, name="shmem_addr_accessible")
            use, intrinsic :: iso_c_binding
            type(c_ptr),         value :: addr
            integer(kind=c_int), value :: pe
        end function shmem_addr_accessible

        integer(kind=c_int) function shmem_my_pe() bind(C, name="shmem_my_pe")
            use, intrinsic :: iso_c_binding
        end function shmem_my_pe

        integer(kind=c_int) function shmem_n_pes() bind(C, name="shmem_n_pes")
            use, intrinsic :: iso_c_binding
        end function shmem_n_pes

        subroutine shmem_barrier_all() bind(C, name="shmem_barrier_all")
            use, intrinsic :: iso_c_binding
        end subroutine shmem_barrier_all

        ! See function signature at
        ! https://cpe.ext.hpe.com/docs/24.07/mpt/openshmemx/shmem_put.html
        ! (accessed 5 March 2025)
        subroutine c_shmem_put(dest, source, nelems, pe) &
#ifdef ENABLE_SHMEM_PUT_GET_32
            bind(C, name="shmem_put32")
#else
            bind(C, name="shmem_putmem")
#endif
            use, intrinsic :: iso_c_binding
            type(c_ptr),            value :: dest
            type(c_ptr),            value :: source
            integer(kind=c_size_t), value :: nelems
            integer(kind=c_int),    value :: pe
        end subroutine c_shmem_put

        subroutine c_shmem_get(dest, source, nelems, pe) &
#ifdef ENABLE_SHMEM_PUT_GET_32
            bind(C, name="shmem_get32")
#else
            bind(C, name="shmem_getmem")
#endif
            use, intrinsic :: iso_c_binding
            type(c_ptr),            value :: dest
            type(c_ptr),            value :: source
            integer(kind=c_size_t), value :: nelems
            integer(kind=c_int),    value :: pe
        end subroutine c_shmem_get
    end interface

contains

    subroutine shmem_logical_malloc(l_data, n)
        use, intrinsic :: iso_c_binding, only : c_ptr, c_sizeof, c_f_pointer
        logical, pointer, intent(inout) :: l_data(:)
        integer,          intent(in)    :: n
        type(c_ptr)                     :: buf_ptr
        logical                         :: l_byte

        buf_ptr = c_shmem_malloc(c_sizeof(l_byte)*n)
        call c_f_pointer(buf_ptr, l_data, [n])

    end subroutine shmem_logical_malloc

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_free(l_data)
        use, intrinsic :: iso_c_binding, only : c_ptr, c_f_pointer, c_loc
        logical, pointer, intent(inout) :: l_data(:)
        type(c_ptr)                     :: buf_ptr

        buf_ptr = c_loc(l_data)
        call c_shmem_free(buf_ptr)
        call c_f_pointer(buf_ptr, l_data, [0])

    end subroutine shmem_free

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_logical_get(dest, source, nelems, pe)
        use, intrinsic :: iso_c_binding, only : c_int, c_size_t, c_ptr, c_loc
#ifndef ENABLE_SHMEM_PUT_GET_32
        use, intrinsic :: iso_c_binding, only : c_sizeof
#endif
        logical, target,  intent(in) :: dest, source
        integer,          intent(in) :: nelems, pe
        integer(c_int)               :: c_pe
        integer(c_size_t)            :: c_nelems
        type(c_ptr)                  :: source_cptr, dest_cptr
#ifndef ENABLE_SHMEM_PUT_GET_32
        logical                      :: l_byte

        c_nelems = nelems * c_sizeof(l_byte)
#else
        c_nelems = nelems
#endif
        c_pe = pe

        source_cptr = c_loc(source)
        dest_cptr = c_loc(dest)

        call c_shmem_get(dest_cptr, source_cptr, c_nelems, c_pe)

    end subroutine shmem_logical_get

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_logical_put(dest, source, nelems, pe)
        use, intrinsic :: iso_c_binding, only : c_int, c_size_t, c_ptr, c_loc
#ifndef ENABLE_SHMEM_PUT_GET_32
        use, intrinsic :: iso_c_binding, only : c_sizeof
#endif
        logical, target,  intent(in) :: dest, source
        integer,          intent(in) :: nelems, pe
        integer(c_int)               :: c_pe
        integer(c_size_t)            :: c_nelems
        type(c_ptr)                  :: source_cptr, dest_cptr
#ifndef ENABLE_SHMEM_PUT_GET_32
        logical                      :: l_byte

        c_nelems = nelems * c_sizeof(l_byte)
#else
        c_nelems = nelems
#endif
        c_pe     = pe

        source_cptr = c_loc(source)
        dest_cptr = c_loc(dest)

        call c_shmem_put(dest_cptr, source_cptr, c_nelems, c_pe)

    end subroutine shmem_logical_put

end module shmem

module parcel_nearest_shmem_graph
    use mpi_layout
    use mpi_utils
    use mpi_timer, only : start_timer       &
                        , stop_timer        &
                        , register_timer
    use parcel_nearest_graph, only : graph_t
    use iso_c_binding, only : c_ptr, c_f_pointer, c_sizeof, c_loc
    use mpi_environment, only : l_ignore_mpi_finalize
    use shmem
    implicit none

    private

    type, extends(graph_t) :: shmem_graph_t

        private

        ! Logicals used to determine which mergers are executed
        ! Integers above could be reused for this, but this would
        ! make the algorithm less readable
        logical, pointer :: l_leaf(:)
        logical, pointer :: l_available(:)
        logical, pointer :: l_merged(:)    ! indicates parcels merged in first stage

        logical :: l_shmem_allocated = .false.

        ! Mapping of neighbouring ranks between MPI Cartesian topology and OpenSHMEM
        integer :: cart2shmem(8)

    contains

        procedure :: initialise => shmem_graph_initialise
        procedure :: finalise   => shmem_graph_finalise
        procedure :: reset      => shmem_graph_reset
        procedure :: resolve    => shmem_graph_resolve
        procedure :: register_timer => shmem_graph_register_timer

        procedure, private :: put_avail
        procedure, private :: put_leaf
        procedure, private :: put_merged

        procedure, private :: get_avail
        procedure, private :: get_leaf
        procedure, private :: get_merged

        procedure, private :: barrier

        procedure, private :: init_cart2shmem
        procedure, private :: get_pe
        procedure, private :: address_check

    end type

    public :: shmem_graph_t

contains

    subroutine shmem_graph_initialise(this, num, l_subcomm)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: num
        logical,              intent(in)    :: l_subcomm
        integer                             :: n, valid
        integer                             :: max_num

        if (this%l_shmem_allocated) then
            return
        endif

        if (.not. l_mpi_layout_initialised) then
            call mpi_stop("Error: The Cartesian communicator not yet initialised.")
        endif

        this%l_shmem_allocated = .true.

#ifndef NULL_ASSIGNMENT_WORKS
        ! this will be set later in parcel_nearest.f90 with a call to create_comm.
        this%comm%comm = MPI_COMM_NULL
#endif

        ! Ensure we use all MPI ranks because we need to call
        ! shmem_barrier_all
        if (l_subcomm) then
            call mpi_print("Ignoring the request to use a subcommunicator.")
            call mpi_print("We only support a global communicator with OpenSHMEM.")
        endif
        this%l_enabled_subcomm = .false.

        call shmem_init

        !--------------------------------------------------
        ! Check mapping beteen MPI and OpenSHMEM

        call this%init_cart2shmem

        !--------------------------------------------------
        ! The routine 'shmem_malloc' must be called with identical
        ! arguments on all processing elements (PEs). This includes the size
        ! of the allocated memory. Otherwise the behaviour of any
        ! following OpenSHMEM call is undefined. Please also see
        ! https://cpe.ext.hpe.com/docs/24.07/mpt/openshmemx/shmem_malloc.html#shmem-malloc
        ! (last accessed 7 March 2025) for further information.
        max_num = num

        call MPI_Allreduce(MPI_IN_PLACE,    &
                           max_num,         &
                           1,               &
                           MPI_INTEGER,     &
                           MPI_MAX,         &
                           world%comm,      &
                           world%err)

        call shmem_logical_malloc(this%l_available, max_num)

        call shmem_logical_malloc(this%l_leaf, max_num)

        call shmem_logical_malloc(this%l_merged, max_num)

        call this%reset


        !--------------------------------------------------
        ! Initial checks to ensure everything works:

        do n = 0, world%size - 1
            if (n /= world%rank) then
                valid = shmem_pe_accessible(n)
                if (valid /= 1) then
                    print *, "Processing element: ", n, " not accessible from ", world%rank
                    call MPI_Abort(MPI_COMM_WORLD, -1, world%err)
                endif

                call this%address_check(this%l_available, n)
                call this%address_check(this%l_leaf, n)
                call this%address_check(this%l_merged, n)
            endif
       enddo

    end subroutine shmem_graph_initialise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_graph_finalise(this)
        class(shmem_graph_t), intent(inout) :: this

        if (.not. this%l_shmem_allocated) then
            return
        endif

        this%l_shmem_allocated = .false.

        call shmem_barrier_all

        call shmem_free(this%l_available)
        call shmem_free(this%l_leaf)
        call shmem_free(this%l_merged)

        call shmem_finalize

        l_ignore_mpi_finalize = .true.

    end subroutine shmem_graph_finalise

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_graph_reset(this)
        class(shmem_graph_t), intent(inout) :: this

        this%l_merged = .false.
        this%l_leaf = .false.
        this%l_available = .false.

    end subroutine shmem_graph_reset

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_graph_resolve(this, isma, iclo, rclo, n_local_small)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(inout) :: isma(0:)
        integer,              intent(inout) :: iclo(:)
        integer,              intent(inout) :: rclo(:)
        integer,              intent(inout) :: n_local_small
        integer                             :: ic, rc, is, m, j
        logical                             :: l_helper
        logical                             :: l_continue_iteration, l_do_merge(n_local_small)
        logical                             :: l_isolated_dual_link(n_local_small)

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
            call this%barrier

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
            call this%barrier

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
            call this%barrier

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
        call this%barrier

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
    end subroutine shmem_graph_resolve

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine shmem_graph_register_timer(this)
        class(shmem_graph_t), intent(inout) :: this

        call this%register_common_timers(label='SHMEM')

    end subroutine shmem_graph_register_timer

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_avail(this, rank, ic, val)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: rank
        integer,              intent(in)    :: ic
        logical,              intent(in)    :: val
        integer                             :: pe

        if (rank == cart%rank) then
            this%l_available(ic) = val
        else
            call start_timer(this%put_timer)
            pe = this%get_pe(rank)
            call shmem_logical_put(dest=this%l_available(ic), source=val, nelems=1, pe=pe)
            call stop_timer(this%put_timer)
        endif

    end subroutine put_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_leaf(this, rank, ic, val)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: rank
        integer,              intent(in)    :: ic
        logical,              intent(in)    :: val
        integer                             :: pe

        if (rank == cart%rank) then
            this%l_leaf(ic) = val
        else
            call start_timer(this%put_timer)
            pe = this%get_pe(rank)
            call shmem_logical_put(dest=this%l_leaf(ic), source=val, nelems=1, pe=pe)
            call stop_timer(this%put_timer)
        endif

    end subroutine put_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine put_merged(this, rank, ic, val)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: rank
        integer,              intent(in)    :: ic
        logical,              intent(in)    :: val
        integer                             :: pe

        if (rank == cart%rank) then
            this%l_merged(ic) = val
        else
            call start_timer(this%put_timer)
            pe = this%get_pe(rank)
            call shmem_logical_put(dest=this%l_merged(ic), source=val, nelems=1, pe=pe)
            call stop_timer(this%put_timer)
        endif

    end subroutine put_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_avail(this, rank, ic) result(val)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: rank
        integer,              intent(in)    :: ic
        logical                             :: val
        integer                             :: pe

        if (rank == cart%rank) then
            val = this%l_available(ic)
        else
            call start_timer(this%get_timer)
            pe = this%get_pe(rank)
            call shmem_logical_get(dest=val, source=this%l_available(ic), nelems=1, pe=pe)
            call stop_timer(this%get_timer)
        endif

    end function get_avail

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_leaf(this, rank, ic) result(val)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: rank
        integer,              intent(in)    :: ic
        logical                             :: val
        integer                             :: pe

        if (rank == cart%rank) then
            val = this%l_leaf(ic)
        else
            call start_timer(this%get_timer)
            pe = this%get_pe(rank)
            call shmem_logical_get(dest=val, source=this%l_leaf(ic), nelems=1, pe=pe)
            call stop_timer(this%get_timer)
        endif

    end function get_leaf

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_merged(this, rank, ic) result(val)
        class(shmem_graph_t), intent(inout) :: this
        integer,              intent(in)    :: rank
        integer,              intent(in)    :: ic
        logical                             :: val
        integer                             :: pe

        if (rank == cart%rank) then
            val = this%l_merged(ic)
        else
            call start_timer(this%get_timer)
            pe = this%get_pe(rank)
            call shmem_logical_get(dest=val, source=this%l_merged(ic), nelems=1, pe=pe)
            call stop_timer(this%get_timer)
        endif

    end function get_merged

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine barrier(this)
        class(shmem_graph_t), intent(inout) :: this

        call start_timer(this%sync_timer)

        call shmem_barrier_all

        call stop_timer(this%sync_timer)

    end subroutine barrier

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_cart2shmem(this)
        class(shmem_graph_t), intent(inout) :: this
        type(MPI_Request)                   :: requests(8)
        type(MPI_Status)                    :: statuses(8)
        integer                             :: n
        integer                             :: me

        ! *this* OpenSHMEM PE (= processing element)
        me = shmem_my_pe()

        ! check if MPI comm world rank == OpenSHMEM rank
        if (world%rank /= me) then
            call mpi_exit_on_error("MPI rank and OpenSHMEM do not agree.")
        endif

        this%cart2shmem = -1

        do n = 1, 8
            call MPI_Isend(me,                      &
                           1,                       &
                           MPI_INTEGER,             &
                           neighbours(n)%rank,      &
                           SEND_NEIGHBOUR_TAG(n),   &
                           cart%comm,               &
                           requests(n),             &
                           cart%err)

            call mpi_check_for_error(cart, &
                "in MPI_Isend of shmem_graph_t::init_cart2shmem.")
        enddo

        do n = 1, 8
            call MPI_Recv(this%cart2shmem(n),       &
                          1,                        &
                          MPI_INTEGER,              &
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
                                "in MPI_Waitall of shmem_graph_t::init_cart2shmem.")

        if (any(this%cart2shmem == -1)) then
            call mpi_exit_on_error("in shmem_graph_t::init_cart2shmem: Not all MPI ranks finished.")
        endif

    end subroutine init_cart2shmem

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    function get_pe(this, rank) result(pe)
        class(shmem_graph_t), intent(in) :: this
        integer,              intent(in) :: rank
        integer                          :: n
        integer                          :: pe

        n = get_neighbour_from_rank(rank)
        pe = this%cart2shmem(n)

    end function get_pe

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine address_check(this, l_array, pe)
        class(shmem_graph_t), intent(in) :: this
        logical, pointer,     intent(in) :: l_array(:)
        integer,              intent(in) :: pe
        integer                          :: valid
        type(c_ptr)                      :: buf_ptr

        buf_ptr = c_loc(l_array)
        valid = shmem_addr_accessible(buf_ptr, pe)
        if (valid /= 1) then
            print *, "Address on PE ", pe, " not accessible from ", world%rank
            call MPI_Abort(MPI_COMM_WORLD, -1, world%err)
        endif

    end subroutine address_check

end module parcel_nearest_shmem_graph
