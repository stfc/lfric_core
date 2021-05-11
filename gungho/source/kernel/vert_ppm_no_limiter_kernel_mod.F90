!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Calculates the coefficients, a0,a1,a2, for 1D subgrid
!>        representation of rho, rho(x) = a0 + a1*x+a2*x**2 with 0<x<1,
!>        in the vertical direction.


!> @details The kernel computes the coefficients a0,a1,a2 where rho is represented
!>          in 1D by the approximation rho(x) = a0+a1*x+a2*x**2
!>          PPM is used to calculate the quadratic subgrid representation of rho.
!>
!>          This kernel is designed to work in the vertical direction only and
!>          takes into account the vertical boundaries.
!>
!>          Note that this kernel only works when rho is a W3 field at lowest order
!>          since it is assumed that ndf_w3 = 1 with stencil_map(1,:) containing
!>          the relevant dofmaps.
module vert_ppm_no_limiter_kernel_mod

use argument_mod,       only : arg_type,          &
                               GH_FIELD, GH_REAL, &
                               GH_READ, GH_WRITE, &
                               CELL_COLUMN
use fs_continuity_mod,  only : W3
use constants_mod,      only : r_def, i_def, l_def
use kernel_mod,         only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vert_ppm_no_limiter_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/             &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: vert_ppm_no_limiter_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vert_ppm_no_limiter_code

contains

!> @brief Compute the subgrid reconstruction coeffiecients for a density field
!! @param[in]  nlayers Number of layers
!! @param[in,out] a0   Coefficient a0
!! @param[in,out] a1   Coefficient a1
!! @param[in,out] a2   Coefficient a2
!! @param[in]  rho     Density
!! @param[in]  ndf_w3  Number of degrees of freedom for W3 per cell
!! @param[in]  undf_w3 Number of unique degrees of freedom for W3
!! @param[in]  map_w3  The dofmap for the cell at the base of the column
subroutine vert_ppm_no_limiter_code( nlayers,                      &
                                     a0,                           &
                                     a1,                           &
                                     a2,                           &
                                     rho,                          &
                                     ndf_w3,                       &
                                     undf_w3,                      &
                                     map_w3 )

  use subgrid_rho_mod, only: second_order_coeffs

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)   :: nlayers
  integer(kind=i_def), intent(in)   :: undf_w3
  real(kind=r_def), intent(inout)   :: a0(undf_w3)
  real(kind=r_def), intent(inout)   :: a1(undf_w3)
  real(kind=r_def), intent(inout)   :: a2(undf_w3)
  real(kind=r_def), intent(in)      :: rho(undf_w3)
  integer(kind=i_def), intent(in)   :: ndf_w3
  integer(kind=i_def), intent(in)   :: map_w3(ndf_w3)

  real(kind=r_def)                  :: coeffs(1:3)
  real(kind=r_def)                  :: rho_local(1:5)
  integer(kind=i_def)               :: k, ii
  logical(kind=l_def)               :: positive, monotone

  ! Apply constant (positive) subgrid reconstruction near the boundaries.
  k = 0
  a0(map_w3(1)+k) = rho(map_w3(1)+k)
  a1(map_w3(1)+k) = 0.0_r_def
  a2(map_w3(1)+k) = 0.0_r_def
  k = nlayers-1
  a0(map_w3(1)+k) = rho(map_w3(1)+k)
  a1(map_w3(1)+k) = 0.0_r_def
  a2(map_w3(1)+k) = 0.0_r_def
  k = 1
  a0(map_w3(1)+k) = rho(map_w3(1)+k)
  a1(map_w3(1)+k) = 0.0_r_def
  a2(map_w3(1)+k) = 0.0_r_def
  k = nlayers-2
  a0(map_w3(1)+k) = rho(map_w3(1)+k)
  a1(map_w3(1)+k) = 0.0_r_def
  a2(map_w3(1)+k) = 0.0_r_def

  positive=.false.
  monotone=.false.

  ! For the layers between 2 and nlayers-3 PPM subgrid reconstruction is applied.
  do k=2,nlayers-3

    do ii=1,5
      rho_local(ii) = rho(map_w3(1)+k+ii-3)
    end do

    call second_order_coeffs(rho_local,coeffs,positive,monotone)
    a0(map_w3(1)+k) = coeffs(1)
    a1(map_w3(1)+k) = coeffs(2)
    a2(map_w3(1)+k) = coeffs(3)

  end do

end subroutine vert_ppm_no_limiter_code

end module vert_ppm_no_limiter_kernel_mod
