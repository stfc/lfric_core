!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the advective update of the physical wind components on the W3
!>        space
!> @details Computes u.grad(v) on scalar points for a specific wind component v.
!>          This is done in the computational space:
!>          u.grad_v = u_bar * (vp-vm) for each wind component (u,v,w)
!>          where u_bar is the wind sampled at the scalar point and
!>          vp & vm are the reconstructed wind on the cell faces
module advective_wind_update_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,   &
                                    GH_FIELD, GH_REAL,     &
                                    GH_READ, GH_WRITE,     &
                                    GH_BASIS, CELL_COLUMN, &
                                    GH_EVALUATOR
use constants_mod,           only : r_def, i_def
use fs_continuity_mod,       only : W2, W3

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: advective_wind_update_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/             &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W2), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W2)  &
       /)
  type(func_type) :: meta_funcs(1) = (/           &
       func_type(W2, GH_BASIS)                    &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: advective_wind_update_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: advective_wind_update_code

contains

!> @param[in] nlayers Number of layers
!> @param[in,out] advective_update u.grad(v) where u is the wind and v is the wind_reconstruction
!> @param[in] wind_reconstruction Reconstructed values of the wind component
!> @param[in] wind The advecting wind field
!> @param[in] ndf_w3 Number of degrees of freedom per cell for the update field
!> @param[in] undf_w3 Number of unique degrees of freedom for the update field
!> @param[in] map_w3 Dofmap for the cell at the base of the column for the update field
!> @param[in] ndf_w2 Number of degrees of freedom per cell for the reconstruction field
!> @param[in] undf_w2 Number of unique degrees of freedom for the reconstruction field
!> @param[in] map_w2 Dofmap for the cell at the base of the column for the reconstruction field
!> @param[in] basis_w2 Basis functions for the vector field evaluated at the
!>                     scalar nodal points
subroutine advective_wind_update_code(nlayers,                 &
                                      advective_update,        &
                                      wind_reconstruction,     &
                                      wind,                    &
                                      ndf_w3, undf_w3, map_w3, &
                                      ndf_w2, undf_w2, map_w2, &
                                      basis_w2                 &
                                      )

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in)    :: nlayers
  integer(kind=i_def),                     intent(in)    :: ndf_w3, undf_w3, &
                                                            ndf_w2, undf_w2
  integer(kind=i_def), dimension(ndf_w3),  intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_w2),  intent(in)    :: map_w2
  real(kind=r_def), dimension(undf_w3),    intent(inout) :: advective_update
  real(kind=r_def), dimension(undf_w2),    intent(in)    :: wind_reconstruction
  real(kind=r_def), dimension(undf_w2),    intent(in)    :: wind

  real(kind=r_def), dimension(3,ndf_w2,ndf_w3), intent(in)  :: basis_w2

  integer(kind=i_def)            :: k, df, df3
  real(kind=r_def), dimension(3) :: grad, uvw

  do k = 0, nlayers-1
    grad(1) = wind_reconstruction(map_w2(3)+k) - wind_reconstruction(map_w2(1)+k)
    grad(2) = wind_reconstruction(map_w2(4)+k) - wind_reconstruction(map_w2(2)+k)
    grad(3) = wind_reconstruction(map_w2(6)+k) - wind_reconstruction(map_w2(5)+k)

    do df3 = 1,ndf_w3
      uvw(:) = 0.0_r_def
      do df = 1,ndf_w2
        uvw(:) = uvw(:) + wind(map_w2(df)+k)*basis_w2(:,df,df3)
      end do
      advective_update(map_w3(df3)+k) = dot_product(uvw, grad)
    end do
  end do

end subroutine advective_wind_update_code

end module advective_wind_update_kernel_mod
