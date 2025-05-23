! =============================================================================
! This module provides the base class for parcel containers.
! =============================================================================
module parcel_container
    use datatypes, only : int64
    use options, only : verbose
    use parameters, only : extent, extenti, center, lower, upper
    use mpi_environment
    use mpi_collectives, only : mpi_blocking_reduce
    use mpi_utils, only : mpi_exit_on_error
    use armanip, only : resize_array
    implicit none

    ! Parcel container type
    type, abstract :: pc_type

        ! number of  parcel attributes
        ! (components are counted individually, e.g. position counts as 3 attributes)
        integer             :: attr_num     ! number of parcel attributes
        integer             :: local_num    ! local number of parcels
        integer(kind=int64) :: total_num    ! global number of parcels (over all MPI ranks)
        integer             :: max_num      ! capacity per attribute, i.e. maximum number of parcels

        ! ---------------------------------------------------------------------
        !   Parcel attributes (common to all types):
        double precision, allocatable, dimension(:, :) :: &
            position,   &
            vorticity,  &
            B               ! B matrix entries; ordering:
                            ! B(:, 1) = B11, B(:, 2) = B12, B(:, 3) = B13
                            ! B(:, 4) = B22, B(:, 5) = B23

        double precision, allocatable, dimension(:) :: &
            volume,     &
            buoyancy

        ! -------------------------------------------------------------------------
        ! buffer indices to access parcel attributes for (de-)serialization;
        ! the buffer indices are set in the extendend types
        integer :: IDX_POS_BEG,         & ! position vector begin
                   IDX_POS_END,         & ! position vector end
                   IDX_VOR_BEG,         & ! vorticity vector begin
                   IDX_VOR_END,         & ! vorticity vector end
                   IDX_SHAPE_BEG,       & ! shape matrix begin
                   IDX_SHAPE_END,       & ! shape matrix end
                   IDX_VOL,             & ! volume
                   IDX_BUO                ! buoyancy

    contains
        ! Base procedures (usually called in derived class procedures):
        procedure :: parcel_base_allocate
        procedure :: parcel_base_deallocate
        procedure :: parcel_base_replace
        procedure :: parcel_base_resize
        procedure :: parcel_base_serialize
        procedure :: parcel_base_deserialize
        ! Basic procedures common to all derived classes:
        procedure :: pack                 => parcel_pack
        procedure :: unpack               => parcel_unpack
        procedure :: delete               => parcel_delete
        ! Pure virtual procedures:
        procedure(parcel_allocate),       deferred :: allocate
        procedure(parcel_deallocate),     deferred :: deallocate
        procedure(parcel_serialize),      deferred :: serialize
        procedure(parcel_deserialize),    deferred :: deserialize
        procedure(parcel_replace),        deferred :: replace
        procedure(parcel_resize),         deferred :: resize
        procedure(parcel_is_small),       deferred :: is_small

    end type pc_type

    interface
        subroutine parcel_allocate(this, num)
            import :: pc_type
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: num
        end subroutine parcel_allocate

        subroutine parcel_deallocate(this)
            import :: pc_type
            class(pc_type), intent(inout) :: this
        end subroutine parcel_deallocate

        subroutine parcel_replace(this, n, m)
            import :: pc_type
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: n, m
        end subroutine

        subroutine parcel_resize(this, new_size)
            import :: pc_type
            class(pc_type), intent(inout) :: this
            integer,        intent(in)    :: new_size
        end subroutine parcel_resize

        subroutine parcel_serialize(this, n, buffer)
            import :: pc_type
            class(pc_type),   intent(in)  :: this
            integer,          intent(in)  :: n
            double precision, intent(out) :: buffer(this%attr_num)
        end subroutine parcel_serialize

        subroutine parcel_deserialize(this, n, buffer)
            import :: pc_type
            class(pc_type),   intent(inout) :: this
            integer,          intent(in)    :: n
            double precision, intent(in)    :: buffer(this%attr_num)
        end subroutine parcel_deserialize

        logical pure function parcel_is_small(this, n)
            import :: pc_type
            class(pc_type), intent(in) :: this
            integer,        intent(in) :: n
        end function parcel_is_small

    end interface

contains

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Allocate parcel memory
    ! ATTENTION: Extended types must allocate additional parcel attributes
    !            in their own routine.
    ! @param[in] num number of parcels
    ! @param[in] n_pos number of spatial dimensions
    ! @param[in] n_vec number of dimensions of vector attributes
    ! @param[in] n_shape number of B matrix elements
    ! @param[in] n_strain number of strain elements
    subroutine parcel_base_allocate(this, num, n_pos, n_vec, n_shape, n_strain)
        class(pc_type), intent(inout) :: this
        integer,        intent(in)    :: num
        integer,        intent(in)    :: n_pos, n_vec, n_shape, n_strain

        allocate(this%position(n_pos, num))
        allocate(this%vorticity(n_vec, num))
        allocate(this%B(n_shape, num))
        allocate(this%volume(num))
        allocate(this%buoyancy(num))

        this%max_num = num

    end subroutine parcel_base_allocate

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Deallocate parcel memory
    ! ATTENTION: Extended types must deallocate additional parcel attributes
    !            in their own routine.
    subroutine parcel_base_deallocate(this)
        class(pc_type), intent(inout) :: this

        if (.not. allocated(this%position)) then
            return
        endif

        this%local_num = 0
        this%total_num = 0
        this%max_num   = 0

        deallocate(this%position)
        deallocate(this%vorticity)
        deallocate(this%B)
        deallocate(this%volume)
        deallocate(this%buoyancy)

    end subroutine parcel_base_deallocate

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Overwrite parcel n with parcel m. This subroutine only replaces the
    ! common types.
    ! ATTENTION: Extended types must replace additional parcel attributes
    !            in their own routine.
    ! @param[in] n index of parcel to be replaced
    ! @param[in] m index of parcel used to replace parcel at index n
    ! @pre n and m must be valid parcel indices
    subroutine parcel_base_replace(this, n, m)
        class(pc_type), intent(inout) :: this
        integer,        intent(in)    :: n, m

        this%position(:, n)  = this%position(:, m)
        this%vorticity(:, n) = this%vorticity(:, m)
        this%volume(n)       = this%volume(m)
        this%buoyancy(n)     = this%buoyancy(m)
        this%B(:, n)         = this%B(:, m)

    end subroutine parcel_base_replace

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Resize the parcel container
    ! ATTENTION: Extended types must resize additional parcel attributes
    !            in their own routine.
    ! @param[in] new_size is the new size of each attribute
    subroutine parcel_base_resize(this, new_size)
        class(pc_type), intent(inout) :: this
        integer,        intent(in)    :: new_size

        if (new_size < this%local_num) then
            call mpi_exit_on_error(&
                "in parcel_container::parcel_base_resize: losing parcels when resizing.")
        endif

        this%max_num = new_size

        call resize_array(this%position, new_size, this%local_num)

        call resize_array(this%vorticity, new_size, this%local_num)
        call resize_array(this%B, new_size, this%local_num)
        call resize_array(this%volume, new_size, this%local_num)
        call resize_array(this%buoyancy, new_size, this%local_num)

    end subroutine parcel_base_resize

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Serialize all parcel attributes into a single buffer
    subroutine parcel_base_serialize(this, n, buffer)
        class(pc_type),   intent(in)  :: this
        integer,          intent(in)  :: n
        double precision, intent(out) :: buffer(this%attr_num)

        buffer(this%IDX_POS_BEG:this%IDX_POS_END)       = this%position(:, n)
        buffer(this%IDX_VOR_BEG:this%IDX_VOR_END)       = this%vorticity(:, n)
        buffer(this%IDX_SHAPE_BEG:this%IDX_SHAPE_END)   = this%B(:, n)
        buffer(this%IDX_VOL)                            = this%volume(n)
        buffer(this%IDX_BUO)                            = this%buoyancy(n)

    end subroutine parcel_base_serialize

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Deserialize all parcel attributes from a single buffer
    subroutine parcel_base_deserialize(this, n, buffer)
        class(pc_type),   intent(inout) :: this
        integer,          intent(in)    :: n
        double precision, intent(in)    :: buffer(this%attr_num)

        this%position(:, n)  = buffer(this%IDX_POS_BEG:this%IDX_POS_END)
        this%vorticity(:, n) = buffer(this%IDX_VOR_BEG:this%IDX_VOR_END)
        this%B(:, n)         = buffer(this%IDX_SHAPE_BEG:this%IDX_SHAPE_END)
        this%volume(n)       = buffer(this%IDX_VOL)
        this%buoyancy(n)     = buffer(this%IDX_BUO)

    end subroutine parcel_base_deserialize

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine parcel_pack(this, pid, num, buffer)
        class(pc_type),   intent(in)  :: this
        integer,          intent(in)  :: pid(:)
        integer,          intent(in)  :: num
        double precision, intent(out) :: buffer(:)
        integer                       :: n, i, j

        do n = 1, num
            i = 1 + (n-1) * this%attr_num
            j = n * this%attr_num
            call this%serialize(pid(n), buffer(i:j))
        enddo
    end subroutine parcel_pack

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine parcel_unpack(this, num, buffer)
        class(pc_type),   intent(inout) :: this
        integer,          intent(in)    :: num
        double precision, intent(in)    :: buffer(:)
        integer                         :: n, i, j

        do n = 1, num
            i = 1 + (n-1) * this%attr_num
            j = n * this%attr_num
            call this%deserialize(this%local_num + n, buffer(i:j))
        enddo

        this%local_num = this%local_num + num

    end subroutine parcel_unpack

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! This algorithm replaces invalid parcels with valid parcels
    ! from the end of the container
    ! @param[in] pid are the parcel indices of the parcels to be deleted
    ! @param[in] n_del is the array size of pid
    ! @pre
    !   - pid must be sorted in ascending order
    !   - pid must be contiguously filled
    !   The above preconditions must be fulfilled so that the
    !   parcel pack algorithm works correctly.
    subroutine parcel_delete(this, pid, n_del)
        class(pc_type), intent(inout) :: this
        integer,        intent(in)    :: pid(0:)
        integer,        intent(in)    :: n_del
        integer                       :: k, l, m

        ! l points always to the last valid parcel
        l = this%local_num

        ! k points always to last invalid parcel in pid
        k = n_del

        ! find last parcel which is not invalid
        do while ((k > 0) .and. (l == pid(k)))
            l = l - 1
            k = k - 1
        enddo

        if (l == -1) then
            call mpi_exit_on_error(&
                "in parcel_container::parcel_delete: more than all parcels are invalid.")
        endif

        ! replace invalid parcels with the last valid parcel
        m = 1

        do while (m <= k)
            ! invalid parcel; overwrite *pid(m)* with last valid parcel *l*
            call this%replace(pid(m), l)

            l = l - 1

            ! find next valid last parcel
            do while ((k > 0) .and. (l == pid(k)))
                l = l - 1
                k = k - 1
            enddo

            ! next invalid
            m = m + 1
        enddo

        ! update number of valid parcels
        this%local_num = this%local_num - n_del

    end subroutine parcel_delete

end module parcel_container
