! =============================================================================
! This module contains global options that can be set at runtime by the user.
! =============================================================================
module options
    implicit none
    !
    ! global options
    !

    ! print more info if true
    logical :: verbose = .false.


    !
    ! output options
    !
    type info
        double precision                  :: parcel_freq        = 1.0d0
        logical                           :: overwrite          = .false.
        character(len=32), dimension(128) :: parcel_list        = ''
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
    end type parcel_type

    type(parcel_type) :: parcel

end module options
