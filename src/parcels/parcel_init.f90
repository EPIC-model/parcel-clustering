! =============================================================================
!               This module initializes parcel default values.
! =============================================================================
module parcel_init
    use options, only : parcel
    use constants, only : zero, two, one, f12, f13, f23, f14
    use parcels_mod, only : parcels
    use parcel_split_mod, only : parcel_split
    use parameters, only : dx, vcell, ncell,            &
                           extent, lower, nx, ny, nz,   &
                           max_num_parcels
    use omp_lib
    use mpi_environment
    use mpi_layout, only : box
    use mpi_utils, only : mpi_print, mpi_exit_on_error
    use mpi_collectives, only : mpi_blocking_reduce
    implicit none

    private :: init_refine

contains

    ! Allocate parcel container and sets values for parcel attributes
    ! to their default values.
    subroutine parcel_default
        double precision             :: lam, l23
        integer                      :: n

        call parcels%allocate(max_num_parcels)

        ! set the number of parcels (see parcels.f90)
        ! we use "n_per_cell" parcels per grid cell
        parcels%local_num = parcel%n_per_cell * box%ncell

        if (parcels%local_num > max_num_parcels) then
            print *, "Number of parcels exceeds limit of", &
                        max_num_parcels, ". Exiting."
            call mpi_exit_on_error
        endif

        parcels%total_num = parcels%local_num
        if (world%size > 1) then
            call mpi_blocking_reduce(parcels%total_num, MPI_SUM, world)
        endif

        call init_regular_positions

        ! initialize the volume of each parcel
        !$omp parallel default(shared)
        !$omp do private(n)
        do n = 1, parcels%local_num
            parcels%volume(n) = vcell / dble(parcel%n_per_cell)
        enddo
        !$omp end do
        !$omp end parallel

        ! aspect ratio: lam = a / c
        lam = maxval((/dx(1) / dx(2), dx(2) / dx(1),   &
                        dx(1) / dx(3), dx(3) / dx(1),   &
                        dx(2) / dx(3), dx(3) / dx(2)/))

        !$omp parallel default(shared)
        !$omp do private(n, l23)
        do n = 1, parcels%local_num
            ! set all to zero
            parcels%B(:, n) = zero

            l23 = (lam * parcels%get_abc(parcels%volume(n))) ** f23

            ! B11
            parcels%B(1, n) = l23

            ! B22
            parcels%B(4, n) = l23
        enddo
        !$omp end do
        !$omp end parallel

        call init_refine(lam)

        !$omp parallel default(shared)
        !$omp do private(n)
        do n = 1, parcels%local_num
            parcels%vorticity(:, n) = zero
            parcels%buoyancy(n) = zero
        enddo
        !$omp end do
        !$omp end parallel

    end subroutine parcel_default

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Position parcels regularly in the domain.
    subroutine init_regular_positions
        integer          :: ix, i, iz, j, iy, k, l, n_per_dim
        double precision :: im, corner(3)

        ! number of parcels per dimension
        n_per_dim = int(dble(parcel%n_per_cell) ** f13)
        if (n_per_dim ** 3 .ne. parcel%n_per_cell) then
            if (world%rank == world%root) then
                print *, "Number of parcels per cell (", &
                            parcel%n_per_cell, ") not a cubic."
            endif
            call mpi_exit_on_error
        endif

        im = one / dble(n_per_dim)

        l = 1
        do iz = 0, nz-1
            do iy = box%lo(2), box%hi(2)
                do ix = box%lo(1), box%hi(1)
                    corner = lower + dble((/ix, iy, iz/)) * dx
                    do k = 1, n_per_dim
                        do j = 1, n_per_dim
                            do i = 1, n_per_dim
                                parcels%position(1, l) = corner(1) + dx(1) * (dble(i) - f12) * im
                                parcels%position(2, l) = corner(2) + dx(2) * (dble(j) - f12) * im
                                parcels%position(3, l) = corner(3) + dx(3) * (dble(k) - f12) * im
                                l = l + 1
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo

        if (.not. parcels%local_num == l - 1) then
            call mpi_exit_on_error("Number of parcels disagree!")
        endif
    end subroutine init_regular_positions

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_refine(lam)
        double precision, intent(inout) :: lam
        double precision                :: evals(3) ! = (a2, b2, c2)

        ! do refining by splitting
        do while (lam >= parcel%lambda_max)
            call parcel_split
            evals = parcels%get_eigenvalues(1)
            lam = sqrt(evals(1) / evals(3))
        end do
    end subroutine init_refine

end module parcel_init
