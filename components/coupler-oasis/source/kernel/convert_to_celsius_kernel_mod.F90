!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Kernel which converts from Kelvin to Celsius
!> @brief but only for temperatures > 1K

module convert_to_celsius_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type,            &
                                    GH_REAL,             &
                                    GH_FIELD, GH_SCALAR, &
                                    GH_WRITE, GH_READ,   &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def

use jules_control_init_mod,        only : n_sea_ice_tile
use lfric_atm_conversions_mod,     only : zero_degrees_celsius

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: convert_to_celsius_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: convert_to_celsius_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: convert_to_celsius_code
contains

!> @brief Kernel which converts from Kelvin to Celsius
!> @brief but only for temperatures > 1K
!! @param[in] nlayers Number of layers
!! @param[in,out] output_field Output field to write to
!! @param[in]     input_field  Input field to read from
!! @param[in]     scalar_a     Scalar field to use in the subtraction
!! @param[in] ndf Number of degrees of freedom per cell for the function space
!! @param[in] undf Number of unique degrees of freedom  for the function space
!! @param[in] map Dofmap for the cell at the base of the column for the function space

subroutine convert_to_celsius_code(nlayers,                     &
                      output_field, input_field,                &
                      ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf
  integer(kind=i_def), intent(in) :: undf
  integer(kind=i_def), dimension(ndf),  intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: output_field
  real(kind=r_def), dimension(undf), intent(in)    :: input_field

  ! Internal variables
  integer(kind=i_def) :: i, df

  do i = 0, n_sea_ice_tile - 1
    do df = 1,ndf
      if (input_field(map(df) + i) > 1.0_r_def) then
        output_field(map(df) + i) = input_field(map(df) + i) -                 &
                                                           zero_degrees_celsius
      else
        output_field(map(df) + i) = 0.0_r_def
      end if
    end do
  end do

end subroutine convert_to_celsius_code

end module convert_to_celsius_kernel_mod
