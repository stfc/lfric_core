!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @details Kernel to initialise IO_Dev global data becase on coordinate value:
!>          this kernel is adapted from nodal_coordinates_kernel, and initialises
!>          each data point to the product of its X, Y and Z coordinates
module io_dev_init_kernel_mod
use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,      &
                                    GH_FIELD, GH_REAL,        &
                                    GH_READ, GH_INC,          &
                                    ANY_SPACE_9, ANY_SPACE_1, &
                                    GH_BASIS, CELL_COLUMN,    &
                                    GH_EVALUATOR
use constants_mod,           only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: io_dev_init_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                       &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  ANY_SPACE_1), &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9)  &
       /)
  type(func_type) :: meta_funcs(1) = (/                     &
       func_type(ANY_SPACE_9, GH_BASIS)                     &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: io_dev_init_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: io_dev_init_code
contains

!> @brief   Initialise IO_Dev data to the product of it's X,Y and Z coordianates
!>          for any function space
!> @param[in]  nlayers   Number of layers
!> @param[in,out] nodal_xyz Nodal coordinate product passed to IO_Dev field
!> @param[in]  chi1      Coordinates in the first direction
!> @param[in]  chi2      Coordinates in the second direction
!> @param[in]  chi3      Coordinates in the third direction
!> @param[in]  ndf_x     Number of degrees of freedom per cell for the output field
!> @param[in]  undf_x    Number of unique degrees of freedom for the output field
!> @param[in]  map_x     Dofmap for the cell at the base of the column for the output field
!> @param[in]  ndf_chi   Number of degrees of freedom per cell for the input field
!> @param[in]  undf_chi  Number of unique degrees of freedom for the input field
!> @param[in]  map_chi   Dofmap for the cell at the base of the column for the input field
!> @param[in]  basis_chi Basis functions of the chi function space evaluated at the
!>                       nodal points of the x function space
subroutine io_dev_init_code( nlayers,                    &
                             nodal_xyz,                  &
                             chi1, chi2, chi3,           &
                             ndf_x, undf_x, map_x,       &
                             ndf_chi, undf_chi, map_chi, &
                             basis_chi )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_x, ndf_chi, undf_x, undf_chi
  integer(kind=i_def), dimension(ndf_x),   intent(in) :: map_x
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  real(kind=r_def), dimension(undf_x),        intent(inout) :: nodal_xyz
  real(kind=r_def), dimension(undf_chi),      intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(1,ndf_chi,ndf_x), intent(in)  :: basis_chi

  ! Internal variables
  integer(kind=i_def) :: df_x, df_chi, k
  real(kind=r_def)    :: xyz(3)

  do k = 0, nlayers-1
    do df_x = 1,ndf_x
      xyz(:) = 0.0_r_def
      do df_chi = 1, ndf_chi
        xyz(1) = xyz(1) + chi1(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(2) = xyz(2) + chi2(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
        xyz(3) = xyz(3) + chi3(map_chi(df_chi)+k)*basis_chi(1,df_chi,df_x)
      end do
      nodal_xyz(map_x(df_x)+k) = xyz(1)*xyz(2)*xyz(3)
    end do
  end do

end subroutine io_dev_init_code

end module io_dev_init_kernel_mod