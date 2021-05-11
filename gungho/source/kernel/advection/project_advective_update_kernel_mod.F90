!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Project an advective update from 3*Ws fields into a W2 field, where
!>        Ws is a scalar function space
!> @details Given the components of a vector stored as 3 W3 fields
!>          compute the Galerkin projection of them into a W2 field.
!>          The use case for this is when the W3 fields contains the components
!>          of the advective update [u.grad(u), u.grad(v), u.grad(w)].
module project_advective_update_kernel_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, func_type,       &
                                    GH_FIELD, GH_REAL,         &
                                    GH_READ, GH_INC,           &
                                    ANY_SPACE_3, ANY_SPACE_9,  &
                                    ANY_DISCONTINUOUS_SPACE_3, &
                                    GH_BASIS, GH_DIFF_BASIS,   &
                                    CELL_COLUMN, GH_QUADRATURE_XYoZ
use constants_mod,           only : r_def, i_def
use fs_continuity_mod,       only : W2

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: project_advective_update_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                    &
       arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_3),              &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
       /)
  type(func_type) :: meta_funcs(3) = (/                                  &
       func_type(W2,          GH_BASIS),                                 &
       func_type(ANY_SPACE_3, GH_BASIS),                                 &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: project_advective_update_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: project_advective_update_code

contains

!> @brief Project an advective update from 3*Ws fields into a W2 field
!! @param[in] nlayers Number of layers
!! @param[in,out] adv_update Projected update into W2
!! @param[in] u_grad_u First component to project: u.grad(u)
!! @param[in] u_grad_v Second component to project: u.grad(v)
!! @param[in] u_grad_w Third component to project: u.grad(w)
!! @param[in] chi1 1st (spherical) coordinate field in Wchi
!! @param[in] chi2 2nd (spherical) coordinate field in Wchi
!! @param[in] chi3 3rd (spherical) coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels.
!! @param[in] ndf_w2 Number of degrees of freedom per cell for vector space
!! @param[in] undf_w2 Number of unique degrees of freedom for vector space
!! @param[in] map_w2 Dofmap for the cell at the base of the column for vector space
!! @param[in] basis_w2 Basis functions for the vector space at quadrature points
!! @param[in] ndf_ws Number of degrees of freedom per cell for scalar space
!! @param[in] undf_ws Number of unique degrees of freedom for scalar space
!! @param[in] map_ws Dofmap for the cell at the base of the column for scalar space
!! @param[in] basis_ws Basis functions for the scalar space at quadrature points
!! @param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!! @param[in] undf_wx Number of unique degrees of freedom for the coordinate space
!! @param[in] map_wx Dofmap for the cell at the base of the column for the coordinate space
!! @param[in] basis_wx Basis functions for the coordinate space at quadrature points
!! @param[in] diff_basis_wx Differential basis functions for the coordinate space at quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Weights of horizontal quadrature points
!! @param[in] wqp_v Weights of vertical quadrature points
subroutine project_advective_update_code(nlayers,                                &
                                         adv_update,                             &
                                         u_grad_u, u_grad_v, u_grad_w,           &
                                         chi1, chi2, chi3,                       &
                                         panel_id,                               &
                                         ndf_w2, undf_w2, map_w2, basis_w2,      &
                                         ndf_ws, undf_ws, map_ws, basis_ws,      &
                                         ndf_wx, undf_wx, map_wx,                &
                                         basis_wx, diff_basis_wx,                &
                                         ndf_pid, undf_pid, map_pid,             &
                                         nqp_h, nqp_v, wqp_h, wqp_v)

  use coordinate_jacobian_mod, only: pointwise_coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def),                     intent(in) :: ndf_w2, ndf_ws, ndf_wx, ndf_pid
  integer(kind=i_def),                     intent(in) :: undf_w2, undf_ws, undf_wx, undf_pid
  integer(kind=i_def), dimension(ndf_ws),  intent(in) :: map_ws
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(3,ndf_w2,nqp_h,nqp_v), intent(in) :: basis_w2
  real(kind=r_def), dimension(1,ndf_ws,nqp_h,nqp_v), intent(in) :: basis_ws
  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v), intent(in) :: basis_wx
  real(kind=r_def), dimension(3,ndf_wx,nqp_h,nqp_v), intent(in) :: diff_basis_wx

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: adv_update
  real(kind=r_def), dimension(undf_ws),  intent(in)    :: u_grad_u, u_grad_v, u_grad_w
  real(kind=r_def), dimension(undf_wx),  intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  real(kind=r_def), intent(in)  :: wqp_h(nqp_h)
  real(kind=r_def), intent(in)  :: wqp_v(nqp_v)

  ! Internal variables
  integer(kind=i_def)                 :: df, k, qp_h, qp_v
  real(kind=r_def), dimension(ndf_wx) :: chi1_e, chi2_e, chi3_e
  real(kind=r_def), dimension(3,3)    :: jac
  real(kind=r_def)                    :: detj
  real(kind=r_def), dimension(3)      :: a3d, v

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1, ndf_wx
      chi1_e(df) = chi1(map_wx(df) + k)
      chi2_e(df) = chi2(map_wx(df) + k)
      chi3_e(df) = chi3(map_wx(df) + k)
    end do

    do qp_v = 1,nqp_v
      do qp_h = 1,nqp_h
        a3d = 0.0_r_def
        do df = 1,ndf_ws
          a3d = a3d + (/ u_grad_u(map_ws(df)+k),  &
                         u_grad_v(map_ws(df)+k),  &
                         u_grad_w(map_ws(df)+k) /)*basis_ws(1,df,qp_h,qp_v)
        end do

        call pointwise_coordinate_jacobian(ndf_wx, chi1_e, chi2_e, chi3_e,  &
                                           ipanel, basis_wx(:,:,qp_h,qp_v), &
                                           diff_basis_wx(:,:,qp_h,qp_v),    &
                                           jac, detj)

        do df = 1,ndf_w2
          ! Advective_update = jac*v . a3d
          v = matmul(jac,basis_w2(:,df,qp_h,qp_v))

          adv_update(map_w2(df)+k) = adv_update(map_w2(df)+k) &
                                   + wqp_h(qp_h)*wqp_v(qp_v)*dot_product(v,a3d)
        end do
      end do
    end do
  end do

end subroutine project_advective_update_code

end module project_advective_update_kernel_mod
