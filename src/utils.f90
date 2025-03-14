module utils
    use constants, only : zero, f12, f23, one, two, pi, twopi
    use parameters, only : lower, vmin, dx, nz, center
    use mpi_timer
    use parcels_mod, only : parcels
    use parcel_merging, only : merge_timer
    use parcel_nearest, only : find_nearest_timer, build_graphs_timer, tree
    use mpi_environment
    use mpi_layout
    use options, only : parcel
    use mpi_utils, only : mpi_exit_on_error, mpi_stop
    implicit none

    private

    integer, allocatable :: seed(:)

    public :: register_all_timers   &
            , init_rng              &
            , setup_parcels

contains

    subroutine register_all_timers

        ! We have 3 + 6 (in trees) timers
        allocate(timings(9))

        call register_timer('parcel merge (total)', merge_timer)
        call register_timer('find nearest', find_nearest_timer)
        call register_timer('build graphs', build_graphs_timer)
    end subroutine register_all_timers

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine init_rng(user_seed)
        integer, intent(in) :: user_seed
        double precision    :: rn
        integer             :: n

        ! 27 Dec 2024
        ! https://stackoverflow.com/a/78623316
        call random_seed(size = n)
        allocate(seed(1:n))
        seed = user_seed + world%rank
        call random_seed(put=seed)

        ! discard first 100 random numbers
        do n = 1, 100
            call random_number(rn)
        enddo

    end subroutine init_rng

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    ! Box-Muller transform
    ! (5 April 2024, https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
    subroutine random_normal(u1, u2)
        double precision, intent(inout) :: u1, u2
        double precision                :: z

        z = sqrt(-two * log(u1))
        u1 = z * cos(twopi * u2)
        u2 = z * sin(twopi * u2)

    end subroutine random_normal

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    ! pick point uniformly on a sphere:
    ! 5 April 2024, https://stats.stackexchange.com/a/7984
    subroutine random_angles(theta, phi)
        double precision, intent(out) :: theta, phi
        double precision              :: u(4)

        ! get 4 uniform numbers in [0, 1[:
        call random_number(u)

        ! transform to standard normal:
        call random_normal(u(1), u(2))
        call random_normal(u(3), u(4))

        ! normalise (note: we do not need u(4) later on):
        u(4) = u(1) ** 2 + u(2) ** 2 + u(3) ** 2
        u(1:3) = u(1:3) / u(4)

        ! azimuthal angle, [0, 2pi[
        theta = atan2(u(2), u(1))

        ! polar angle, [0, pi[
        u(3) = max(-1.0d0, min(u(3), 1.0d0))
        phi = acos(u(3))

    end subroutine random_angles

    !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

    subroutine setup_parcels(ratio, l_shuffle, l_variable_nppc)
        double precision, intent(in) :: ratio
        logical,          intent(in) :: l_shuffle, l_variable_nppc
        double precision             :: rn(10), lam, lam2, abc, a2, b2, c2, theta, phi
        double precision             :: st, ct, sp, cp, corner(3), x, y, z
        integer                      :: ix, iy, iz, m, l, npp, n_per_cell

        npp = parcel%n_per_cell
        n_per_cell = max(10, parcel%n_per_cell - 10)

        if ((ratio < 0.0d0) .or. (ratio > 1.0d0)) then
            call mpi_stop("Fraction of small parcels must be in [0, 1].")
        endif

        l = 1
        do iz = 0, nz-1
            do iy = box%lo(2), box%hi(2)
                do ix = box%lo(1), box%hi(1)
                    if (l_variable_nppc) then
                        call random_number(lam)
                        npp = int(lam * n_per_cell) + 10    ! ensure at least 10 parcels per cell
                    endif
                    corner = lower + dble((/ix, iy, iz/)) * dx
                    do m = 1, npp
                        ! rn between 0 and 1
                        call random_number(rn)

                        x = corner(1) + dx(1) * rn(1)
                        y = corner(2) + dx(2) * rn(2)
                        z = corner(3) + dx(3) * rn(3)

                        parcels%position(1, l) = x
                        parcels%position(2, l) = y
                        parcels%position(3, l) = z

                        ! vorticity between -10 and 10: y = 20 * x - 10
                        parcels%vorticity(1, l) = 20.0d0 * rn(4) - 10.d0
                        parcels%vorticity(2, l) = 20.0d0 * rn(5) - 10.d0
                        parcels%vorticity(3, l) = 20.0d0 * rn(6) - 10.d0

                        ! buoyancy between -1 and 1: y = 2 * x - 1
                        parcels%buoyancy(l) = 2.0d0 * rn(7) - 1.d0

                        ! volume between [0, 1] * vmin + (1 - r) * vmin
                        ! where r in [0, 1] is the fraction of small parcels
                        ! we ensure that no parcel is smaller than 1/4 of vmin
                        parcels%volume(l) = max(0.25d0 * vmin, vmin * rn(8) + (1.0d0 - ratio) * vmin)

                        ! lam = a / c in [1, 4]
                        lam = 3.d0 * rn(9) + 1.0d0

                        ! lam2 = a / b
                        lam2 = 3.d0 * rn(10) + 1.0d0

                        abc = 0.75d0 * parcels%volume(l) / pi

                        a2 = (abc * lam * lam2)  ** f23
                        b2 = a2 / lam2 ** 2
                        c2 = a2 / lam ** 2

                        ! get random angles:
                        call random_angles(theta, phi)

                        st = sin(theta)
                        ct = cos(theta)
                        sp = sin(phi)
                        cp = cos(phi)

                        parcels%B(1, l) = a2 * ct ** 2 * sp ** 2 + b2 * st ** 2 + c2 * ct ** 2 * cp ** 2
                        parcels%B(2, l) = a2 * st * ct * sp ** 2 - b2 * st * ct + c2 * st * ct * cp ** 2
                        parcels%B(3, l) = (a2 - c2) * ct * sp * cp
                        parcels%B(4, l) = a2 * st ** 2 * sp ** 2 + b2 * ct ** 2 + c2 * st ** 2 * cp ** 2
                        parcels%B(5, l) = (a2 - c2) * st * sp * cp

                        l = l + 1
                    enddo
                enddo
            enddo
        enddo
        parcels%local_num = l - 1

        if (l_shuffle) then
            call shuffleall
        endif

    end subroutine setup_parcels

    ! performs a Fisher-Yates shuffle (aka Knuth shuffle)
    subroutine shuffleall
        integer          :: shuffle_index, rand_target
        double precision :: tmp_var
        double precision :: tmp_vec(3), tmp_B(5)
        double precision :: random_out

        do shuffle_index = parcels%local_num, 2, -1
            call random_number(random_out)
            rand_target = int(random_out * shuffle_index) + 1

            tmp_vec = parcels%position(:, rand_target)
            parcels%position(:, rand_target) = parcels%position(:, shuffle_index)
            parcels%position(:, shuffle_index) = tmp_vec

            tmp_vec = parcels%vorticity(:, rand_target)
            parcels%vorticity(:, rand_target) = parcels%vorticity(:, shuffle_index)
            parcels%vorticity(:, shuffle_index) = tmp_vec

            tmp_var = parcels%buoyancy(rand_target)
            parcels%buoyancy(rand_target) = parcels%buoyancy(shuffle_index)
            parcels%buoyancy(shuffle_index) = tmp_var

            tmp_var = parcels%volume(rand_target)
            parcels%volume(rand_target) = parcels%volume(shuffle_index)
            parcels%volume(shuffle_index) = tmp_var

            tmp_B = parcels%B(:, rand_target)
            parcels%B(:, rand_target) = parcels%B(:, shuffle_index)
            parcels%B(:, shuffle_index) = tmp_B
        end do
    end subroutine shuffleall

end module utils
