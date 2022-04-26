!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Process ocean and sea ice data coming through the coupler

module process_o2a_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,              &
                                    GH_FIELD, GH_REAL,     &
                                    GH_READ, GH_READWRITE, &
                                    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def
use surface_config_mod,      only : therm_cond_sice,       &
                                    therm_cond_sice_snow
use jules_control_init_mod,        only : n_sea_ice_tile
use lfric_atm_conversions_mod,     only: zero_degrees_celsius

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: process_o2a_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                    &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: process_o2a_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: process_o2a_code
contains

!> @brief Process ocean and sea ice data coming through the coupler
!! @param[in] nlayers Number of layers
!! @param[in,out] sea_ice_fraction     the sea ice fraction field
!! @param[in,out] sea_ice_thickness    the sea ice thickness field
!! @param[in,out] sea_ice_layer_t      the sea ice layer temperature field
!! @param[in] ndf Number of degrees of freedom per cell for the function space
!! @param[in] undf Number of unique degrees of freedom  for the function space
!! @param[in] map Dofmap for the cell at the base of the column for the function space

subroutine process_o2a_code(nlayers,                                 &
                              sea_ice_fraction, sea_ice_thickness,   &
                              sea_ice_layer_t,                       &
                              ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf
  integer(kind=i_def), intent(in) :: undf
  integer(kind=i_def), dimension(ndf),  intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: sea_ice_fraction
  real(kind=r_def), dimension(undf), intent(inout) :: sea_ice_thickness
  real(kind=r_def), dimension(undf), intent(inout) :: sea_ice_layer_t

  ! Parameters
  ! Minimum sea ice thickness allowed through the coupler
  real(kind=r_def), parameter :: hi_min = 0.05_r_def
  ! Minimum sea ice fraction allowed through the coupler
  real(kind=r_def), parameter :: aicenmin = 1.0e-02_r_def
  ! Maximum sea ice layer temperature allowed through the coupler
  real(kind=r_def), parameter :: ti_max = zero_degrees_celsius
  ! Minimum sea ice layer temperature allowed through the coupler
  real(kind=r_def), parameter :: ti_min = 200.0_r_def

  ! Local variables
  integer(kind=i_def) :: i
  real(kind=r_def)    :: total_sea_ice_fraction



  total_sea_ice_fraction = 0.0_r_def

  do i = 0, n_sea_ice_tile - 1

      ! Make sure sea ice is thicker than hi_min
      if ( sea_ice_fraction(map(1) + i) > 0.0_r_def ) then
        if ( sea_ice_thickness(map(1) + i) / sea_ice_fraction(map(1) + i) < hi_min ) then
          sea_ice_fraction(map(1) + i) = 0.0_r_def
        end if
      end if

      if ( sea_ice_fraction(map(1) + i)  < aicenmin ) then
        sea_ice_fraction(map(1) + i) = 0.0_r_def
        sea_ice_thickness(map(1) + i) = 0.0_r_def
        sea_ice_layer_t(map(1) + i) = 0.0_r_def
      end if

      total_sea_ice_fraction = total_sea_ice_fraction + sea_ice_fraction(map(1) + i)

  end do

  do i = 0, n_sea_ice_tile - 1

      ! Reduce category ice fractions if aggregate > 1
      if (total_sea_ice_fraction > 1.0_r_def) then
        sea_ice_fraction(map(1) + i) = sea_ice_fraction(map(1) + i) /          &
                                                         total_sea_ice_fraction
      end if

      ! Undo sea ice fraction scaling done on ocean side of coupler
      if ( sea_ice_fraction(map(1) + i) > 0.0_r_def ) then
        sea_ice_thickness(map(1) + i) = sea_ice_thickness(map(1) + i) /        &
                                                   sea_ice_fraction(map(1) + i)
        sea_ice_layer_t(map(1) + i) = sea_ice_layer_t(map(1) + i) /            &
                                                   sea_ice_fraction(map(1) + i)
      end if

      ! Apply bounds to sea ice layer temperatures
      if ( sea_ice_layer_t(map(1) + i) > ti_max ) then
        sea_ice_layer_t(map(1) + i) = ti_max
      end if
      if ( sea_ice_layer_t(map(1) + i) < ti_min .and.                          &
                                sea_ice_layer_t(map(1) + i) /= 0.0_r_def ) then
        sea_ice_layer_t(map(1) + i) = ti_min
      end if
  end do

end subroutine process_o2a_code

end module process_o2a_kernel_mod
