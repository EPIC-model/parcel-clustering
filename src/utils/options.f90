! =============================================================================
! This module contains global options that can be set at runtime by the user.
! =============================================================================
module options
    use mpi_utils, only : mpi_stop
    implicit none
    !
    ! global options
    !

    ! print more info if true
    logical :: verbose = .false.

    ! if a restarted simulation
    logical :: l_restart = .false.

    ! configuration file
    character(len=512) :: filename = ''

    ! restart file
    character(len=512) :: restart_file = ''

    ! field input file
    character(len=512) :: field_file = ''


    !
    ! output options
    !
    type info
        double precision                  :: field_freq         = 1.0d0
        logical                           :: write_fields       = .true.
        character(len=32), dimension(128) :: field_list         = ''
        double precision                  :: parcel_freq        = 1.0d0
        logical                           :: overwrite          = .false.
        logical                           :: write_parcels      = .true.
        character(len=32), dimension(128) :: parcel_list        = ''
        double precision                  :: parcel_stats_freq  = 1.0d0
        logical                           :: write_parcel_stats = .true.
        double precision                  :: field_stats_freq   = 1.0d0
        logical                           :: write_field_stats  = .true.
        character(len=512)                :: basename           = ''
    end type info

    type(info) :: output

    !
    ! parcel options
    !
    type parcel_type
        double precision :: size_factor      = 1.0d0    ! factor to increase max. number of parcels
        double precision :: grow_factor      = 1.2d0    ! factor to increase the parcel container size
                                                        ! in the parcel splitting routine
        double precision :: shrink_factor    = 0.8d0    ! factor to reduce the parcel container size
        integer          :: n_per_cell       = 8        ! number of parcels per cell (need to be a cube)
        double precision :: lambda_max       = 4.0d0    ! max. ellipse aspect ratio a/b
        double precision :: min_vratio       = 20.0d0   ! minimum ratio of grid cell volume / parcel volume
        integer          :: correction_iters = 2        ! parcel correction iterations
        double precision :: gradient_pref    = 1.8d0    ! prefactor for gradient descent
        double precision :: max_compression  = 0.5d0    ! parameter for gradient descent
                                                        ! (limits the shift in parcel position)
    end type parcel_type

    type(parcel_type) :: parcel

contains
    ! parse configuration file
    ! (see https://cyber.dabamos.de/programming/modernfortran/namelists.html [8 March 2021])
    subroutine read_config_file
        integer :: ios
        integer :: fn = 1
        logical :: exists = .false.

        ! namelist definitions
        namelist /EPIC/ field_file, output, parcel

        ! check whether file exists
        inquire(file=filename, exist=exists)

        if (exists .eqv. .false.) then
            call mpi_stop(&
                'Error: input file "' // trim(filename) // '" does not exist.')
        endif

        ! open and read Namelist file.
        open(action='read', file=filename, iostat=ios, newunit=fn)

        read(nml=EPIC, iostat=ios, unit=fn)

        if (ios /= 0) then
            call mpi_stop('Error: invalid Namelist format.')
        end if

        close(fn)

        ! check whether NetCDF files already exist
        inquire(file=output%basename, exist=exists)

        if (exists) then
            call mpi_stop(&
                'Error: output file "' // trim(output%basename) // '" already exists.')
        endif

    end subroutine read_config_file

end module options
