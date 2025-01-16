! =============================================================================
! This module contains global parameters that stay constant throughout a
! simulation.
! =============================================================================
module parameters
    use options, only : parcel
    use datatypes, only : int64
    use constants
    use netcdf_reader
    use netcdf_utils
    use netcdf_writer
    use mpi_environment
    use mpi_layout, only : box, l_mpi_layout_initialised
    use mpi_utils, only : mpi_print, mpi_stop
    implicit none

    ! mesh spacing
    double precision, protected :: dx(3)

    ! inverse mesh spacing
    double precision, protected :: dxi(3)

    ! grid cell volume
    double precision, protected :: vcell

    ! inverse grid cell volume
    double precision, protected :: vcelli

    ! number of grid cells in each dimension
    integer :: nx, ny, nz

    ! total number of grid cells
    integer(kind=int64), protected :: ncell

    ! inverse of total number of grid cells
    double precision, protected :: ncelli

    ! total number of grid cells in horizontal plane (x, y)
    integer, protected :: nhcell

    ! inverse of total number of grid cells in horizontal plane (x, y)
    double precision, protected :: nhcelli

    ! total number of grid points
    integer, protected :: ngrid

    ! inverse of total number of grid points
    double precision, protected :: ngridi

    ! domain size
    double precision :: extent(3)

    ! inverse domain size
    double precision, protected :: extenti(3)

    ! domain volume
    double precision, protected :: vdomain

    ! inverse domain volume
    double precision, protected :: vdomaini

    ! domain centre
    double precision, protected :: center(3)

    ! domain half widths values
    double precision, protected :: hl(3)

    double precision, protected :: hli(3)

    ! domain origin
    double precision :: lower(3)

    ! domain upper boundary
    double precision, protected :: upper(3)

    ! minimum volume
    double precision, protected :: vmin

    ! upper bound for major semi-axis (used for splitting)
    double precision, protected :: amax

    ! maximum number of allowed parcels
    integer, protected :: max_num_parcels

    ! specifies if zeta is kept zero on a boundary;
    ! this also makes sure that dzeta/dt = 0 on a boundary
    logical, protected :: l_bndry_zeta_zero(2)

contains

    ! Update all parameters according to the
    ! user-defined global options.
    subroutine update_parameters
        double precision    :: msr
        integer(kind=int64) ::  max_size

        if (.not. l_mpi_layout_initialised) then
            call mpi_print("MPI layout is not initialsed!")
            call mpi_stop
        endif

        upper = lower + extent

        extenti = one / extent

        dx = extent / dble((/nx, ny, nz/))
        dxi = one / dx;

        msr = maxval((/dxi(1) * dx(2), dxi(2) * dx(1),   &
                       dxi(1) * dx(3), dxi(3) * dx(1),   &
                       dxi(2) * dx(3), dxi(3) * dx(2)/))

        if (msr > two) then
            if (world%rank == world%root) then
                print *, '**********************************************************************'
                print *, '*                                                                    *'
                print *, '*   Warning: A mesh spacing ratio of more than 2 is not advisable!   *'
                print *, '*                                                                    *'
                print *, '**********************************************************************'
            endif
        endif

        vdomain = product(extent)
        vdomaini = one / vdomain

        vcell = product(dx)
        vcelli = one / vcell

        nhcell = nx * ny
        nhcelli = one / dble(nhcell)

        ncell = nhcell * nz
        ncelli = one / dble(ncell)

        ! due to x periodicity it is only nx
        ngrid = nx * ny * (nz + 1)
        ngridi = one / dble(ngrid)

        ! domain
        center = f12 * (lower + upper)
        hl = extent / two
        hli = one / hl

        vmin = vcell / parcel%min_vratio

        amax = (f34 * fpi) ** f13 * minval(dx)

        max_size = int(box%halo_ncell * parcel%min_vratio * parcel%size_factor, kind=int64)

        if (max_size > huge(max_num_parcels)) then
            if (world%rank == world%root) then
                print *, "Error: Maximum number of parcels larger than integer"
                print *, "       representation. Overflow! You can circumvent this"
                print *, "       issue by using more MPI cores."
            endif
            call mpi_stop
        endif

        max_num_parcels = int(max_size)

    end subroutine update_parameters

    subroutine set_mesh_spacing(ext, nc)
        double precision, intent(in) :: ext(3)
        integer,          intent(in) :: nc(3)
        dx = ext / dble(nc)
    end subroutine set_mesh_spacing


    subroutine set_max_num_parcels(num)
        integer, intent(in) :: num
        max_num_parcels = num
    end subroutine set_max_num_parcels

end module parameters
