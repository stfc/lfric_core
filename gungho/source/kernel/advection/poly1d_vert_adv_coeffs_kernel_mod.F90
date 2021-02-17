!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the coefficients for reconstructing a
!>        1D vertical upwind polynomial representation of a tracer field on the
!>        centres of a cell
!> @details Compute the coefficients to reconstruct a tracer field using a high order
!>          1D polynomial fit to the values of the tracer over a given stencil.
!>          The stencil used for the polynomial is centred on the upwind cell for each edge
!>          A polynomial is used containing all monomials up to the
!>          desired order, i.e. order = 2: 1 + x + x^2
!>          This is exactly fitted over all cells in the stencil
!>          The methodology is inspired by that of Thuburn et.al GMD 2014 for
!>          2D reconstructions
!>          This method is only valid for lowest order elements
module poly1d_vert_adv_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type,  &
                              GH_FIELD, GH_INTEGER, &
                              GH_WRITE, GH_READ,    &
                              ANY_SPACE_1,          &
                              GH_BASIS, CELLS, GH_EVALUATOR
use constants_mod,     only : r_def, i_def, EPS
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_adv_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, Wtheta),                         &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_1),                    &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: poly1d_vert_adv_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_adv_coeffs_code
contains

!>@brief Compute the coefficients needed for a 1D vertical reconstruction
!>       of a tracer field on vertical faces
!>@param[in] nlayers Number of vertical layers
!>@param[out] coeff Array of fields to store the coefficients for the polynomial
!!                  reconstruction
!>@param[in] chi1 1st component of the physical coordinate field
!>@param[in] chi2 2nd component of the physical coordinate field
!>@param[in] chi3 3rd component of the physical coordinate field
!>@param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!>@param[in] undf_wt Total number of degrees of freedom for Wtheta
!>@param[in] map_wt Dofmap of the tracer field
!>@param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!>@param[in] undf_wx Total number of degrees of freedom for the coordinate space
!>@param[in] map_wx Dofmap of the coordinate space
!>@param[in] basis_wx Basis function of the coordinate space evaluated on
!!                    Wtheta nodal points
!>@param[in] global_order Desired polynomial order for advective computations
!>@param[in] nfaces_v Number of vertical faces (used by PSyclone to size coeff
!!                    array)
subroutine poly1d_vert_adv_coeffs_code(nlayers,                   &
                                       coeff,                     &
                                       chi1, chi2, chi3,          &
                                       ndf_wt,                    &
                                       undf_wt,                   &
                                       map_wt,                    &
                                       ndf_wx,                    &
                                       undf_wx,                   &
                                       map_wx,                    &
                                       basis_wx,                  &
                                       global_order,              &
                                       nfaces_v )


  use matrix_invert_mod,         only: matrix_invert
  use base_mesh_config_mod,      only: geometry, &
                                       geometry_spherical
  use finite_element_config_mod, only: spherical_coord_system,     &
                                       spherical_coord_system_xyz
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: global_order, nfaces_v
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, &
                                     ndf_wx, undf_wx

  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx

  real(kind=r_def), dimension(undf_wx), intent(in) :: chi1, chi2, chi3

  real(kind=r_def), dimension(global_order+1, nfaces_v, undf_wt), intent(out) :: coeff

  real(kind=r_def), dimension(1,ndf_wx,ndf_wt), intent(in) :: basis_wx

  ! Local variables
  integer(kind=i_def) :: k, kk, ijk, df, stencil, nmonomial, qp, &
                         m, direction, order, kx, spherical, planar, &
                         boundary_offset, vertical_order
  integer(kind=i_def), allocatable, dimension(:,:)         :: smap
  real(kind=r_def)                                         :: xx, fn, z0, zq
  real(kind=r_def),              dimension(3)              :: x0, xq
  real(kind=r_def), allocatable, dimension(:,:)            :: monomial_matrix, &
                                                              inv_monomial_matrix
  real(kind=r_def)                                         :: basis_at_w3

  ! Ensure that we reduce the order if there are only a few layers
  vertical_order = min(global_order, nlayers-1)

  if ( geometry == geometry_spherical .and. &
       spherical_coord_system == spherical_coord_system_xyz ) then
    spherical = 1.0_r_def
    planar    = 0.0_r_def
  else
    spherical = 0.0_r_def
    planar    = 1.0_r_def
  end if

  ! Compute the offset map for all orders up to order
  allocate( smap(global_order+1,0:global_order) )
  smap(:,:) = 0
  do m = 0,global_order
    do stencil = 1,m+1
      smap(stencil,m) = - m/2 + (stencil-1)
    end do
  end do

  ! Step 1: Build monomials over all cells in advection stencils

  ! Loop over layers
  layer_loop: do k = 0, nlayers-1

   ! Position vector in the centre of this cell
    x0 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = map_wx( df ) + k
      basis_at_w3 = 0.5_r_def*(basis_wx(1,df,1) + basis_wx(1,df,2))
      x0(:) = x0(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_at_w3
    end do
    z0 = sqrt(x0(1)**2 + x0(2)**2 + x0(3)**2)*spherical + x0(3)*planar

    ! Compute the coefficients of each point in the stencil for
    ! each cell when this is the upwind point
    ! Loop over both posivie and negative directions
    direction_loop: do direction = 0,1

      ! Compute local order, this is at most vertical_order but reduces near the
      ! top and bottom boundary
      order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))

      ! For the boundary points reduce to centred linear interpolation
      if ( ( k == 0           .and. direction == 0) .or. &
           ( k == nlayers - 1 .and. direction == 1) ) order = 1

      ! Number of monomials to use (all polynomials up to total degree of order+1)
      nmonomial = (order + 1)
      allocate( monomial_matrix(nmonomial, nmonomial),  &
                inv_monomial_matrix(nmonomial, nmonomial) )

      ! Offset term to make sure we use the right indices at the top
      ! boundary
      boundary_offset = 0
      if ( k == nlayers - 1 .and. direction == 1)  boundary_offset = -1
      monomial_matrix = 0.0_r_def

      ! Loop over all cells in the stencil
      stencil_loop: do stencil = 1, order+1
        kk = k + smap(stencil,order) + direction + boundary_offset
        ! If k + smap(stencil,order/2) == nlayers we need to make sure we use
        ! coordinate evaluated in cell k = nlayers-1 but with zp=1
        if ( kk == nlayers ) then
          kx = nlayers-1
          qp = 2
        else
          kx = kk
          qp = 1
        end if
        ! Evaluate coordinate at theta point of this cell
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = map_wx( df ) + kx
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qp)
        end do
        zq = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)*spherical + xq(3)*planar

        ! Second: Compute the local coordinate of each quadrature point from the
        !         physical coordinate
        xx = (zq - z0)

        ! Third: Compute each needed monomial in terms of the local coordinate
        !        on each quadrature point
        ! Loop over monomials
        do m = 1, nmonomial
          fn = xx**(m-1)
          monomial_matrix(stencil,m) = monomial_matrix(stencil,m) + fn
        end do
      end do stencil_loop

      ! Manipulate the integrals of monomials,
      call matrix_invert(monomial_matrix,inv_monomial_matrix,nmonomial)

      ! Initialise polynomial coeffficients to zero
      coeff(:,direction+1,map_wt(1)+k) = 0.0_r_def

      ! Evaluate the polynomial at Wtheta point which is the origin of our local
      ! coordinate system (x=0),
      ! The process for doing this would be to compute momonmial(j) = x^(j-1)
      ! then compute coefficients as:
      ! C_i = dot_product(monomial, inv_monomial_matrix * delta)
      ! with delta_i = 1, delta_j = 0, j /= i
      !
      ! but since x = 0, we have monomial(:) = (1,0,..,0)
      ! and so the above is simplified to
      ! C_i = inv_monomial_matrix(1,i)
      do stencil = 1, order+1
        coeff(stencil,direction+1,map_wt(1)+k) = inv_monomial_matrix(1,stencil)
      end do
      deallocate( monomial_matrix, inv_monomial_matrix )
    end do direction_loop
  end do layer_loop
  ! Coeffecients are really for W3 points, but we are using a Wtheta field
  ! which contains an extra point so just set those to be zero
  coeff(:,:,map_wt(2) + nlayers-1 ) = 0.0_r_def
  deallocate( smap )

end subroutine poly1d_vert_adv_coeffs_code

end module poly1d_vert_adv_coeffs_kernel_mod
