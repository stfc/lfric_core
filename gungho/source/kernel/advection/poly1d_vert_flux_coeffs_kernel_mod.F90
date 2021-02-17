!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Compute the coefficients for reconstructing a
!>        a  1D vertical upwind polynomial representation of a tracer field on the
!>        faces of a cell
!> @details Compute the coefficients of the flux of a tracer density field using a high order
!>          1D polynomial fit to the integrated tracer values over a given stencil.
!>          The stencil used for the polynomial is centred on the upwind cell for each edge
!>          A symmetric polynomial is used containing all monomials up to the
!>          desired order, i.e. order = 2: 1 + x + x^2
!>          This is exactly fitted over all cells in the stencil
!>          The methodology is inspired by that of Thuburn et.al GMD 2014 for
!>          2D reconstructions
!>          This method is only valid for lowest order elements
module poly1d_vert_flux_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type,  &
                              GH_FIELD, GH_INTEGER, &
                              GH_WRITE, GH_READ,    &
                              ANY_SPACE_1,          &
                              GH_BASIS, CELLS,      &
                              GH_QUADRATURE_XYoZ,   &
                              GH_QUADRATURE_face
use constants_mod,     only : r_def, i_def, EPS
use fs_continuity_mod, only : W3
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly1d_vert_flux_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W3),                             &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_1),                    &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape(2) = (/ GH_QUADRATURE_XYoZ, GH_QUADRATURE_face /)
contains
  procedure, nopass :: poly1d_vert_flux_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly1d_vert_flux_coeffs_code
contains

!>@brief Compute the coefficients needed for a 1D vertical reconstruction
!>       of a tracer field on vertical faces
!>@param[in] nlayers Number of vertical layers
!>@param[out] coeff Array of fields to store the coefficients for the polynomial
!!                  reconstruction
!>@param[in]  mdw3 Mass matrix diagonal for the W3 space, this is used to give
!!                 the cell volume
!>@param[in] chi1 1st component of the physical coordinate field
!>@param[in] chi2 2nd component of the physical coordinate field
!>@param[in] chi3 3rd component of the physical coordinate field
!>@param[in] ndf_w3 Number of degrees of freedom per cell for W3
!>@param[in] undf_w3 Total number of degrees of freedom for W3
!>@param[in] map_w3 Dofmap of the tracer field
!>@param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!>@param[in] undf_wx Total number of degrees of freedom for the coordinate space
!>@param[in] map_wx Dofmap of the coordinate space stencil
!>@param[in] basis_wx Basis function of the coordinate space evaluated on
!!                    quadrature points
!>@param[in] face_basis_wx Basis function of the coordinate space evaluated on
!!                         quadrature points on the vertical faces
!>@param[in] global_order Desired polynomial order for flux computations
!>@param[in] nfaces_re_v Number of vertical faces( used by PSyclone to size
!!                       coeff array)
!>@param[in] nqp_h Number of horizontal quadrature points
!>@param[in] nqp_v Number of vertical quadrature points
!>@param[in] wqp_h Weights of horizontal quadrature points
!>@param[in] wqp_v Weights of vertical quadrature points
!>@param[in] nfaces_qr Number of faces in the quadrature rule
!>@param[in] nqp_f Number of face quadrature points
!>@param[in] wqp_f Weights of face quadrature points
subroutine poly1d_vert_flux_coeffs_code(nlayers,                    &
                                        coeff,                      &
                                        mdw3,                       &
                                        chi1, chi2, chi3,           &
                                        ndf_w3,                     &
                                        undf_w3,                    &
                                        map_w3,                     &
                                        ndf_wx,                     &
                                        undf_wx,                    &
                                        map_wx,                     &
                                        basis_wx,                   &
                                        face_basis_wx,              &
                                        global_order,               &
                                        nfaces_re_v,                &
                                        nqp_h, nqp_v, wqp_h, wqp_v, &
                                        nfaces_qr, nqp_f, wqp_f )


  use matrix_invert_mod,          only: matrix_invert
  use base_mesh_config_mod,       only: geometry, &
                                        geometry_spherical
  use finite_element_config_mod,  only: spherical_coord_system,     &
                                        spherical_coord_system_xyz
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: global_order, nfaces_re_v
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3, &
                                     ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: nqp_v, nqp_h, nqp_f, nfaces_qr

  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx

  real(kind=r_def), dimension(undf_w3), intent(in)  :: mdw3
  real(kind=r_def), dimension(undf_wx), intent(in)  :: chi1, chi2, chi3

  real(kind=r_def), dimension(global_order+1, nfaces_re_v, undf_w3), intent(out) :: coeff

  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v),     intent(in) :: basis_wx
  real(kind=r_def), dimension(1,ndf_wx,nqp_f,nfaces_qr), intent(in) :: face_basis_wx

  real(kind=r_def), dimension(nqp_h),           intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v),           intent(in) ::  wqp_v
  real(kind=r_def), dimension(nqp_f,nfaces_qr), intent(in) ::  wqp_f

  ! Local variables
  integer(kind=i_def) :: k, ijk, df, qv0, qh0, stencil, nmonomial, qp, &
                         m, face, order, boundary_offset, offset, &
                         vertical_order
  integer(kind=i_def), allocatable, dimension(:,:)         :: smap
  real(kind=r_def)                                         :: xx, fn, z0, zq, spherical, planar
  real(kind=r_def),              dimension(3)              :: x0, xq
  real(kind=r_def), allocatable, dimension(:,:)            :: int_monomial, inv_int_monomial
  real(kind=r_def), allocatable, dimension(:)              :: beta, delta, monomial
  real(kind=r_def),              dimension(global_order+1) :: area

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

  ! Avoid compile warning for unused variables
  fn = wqp_h(1)
  fn = wqp_f(1,1)

  ! Index of quadrature point in the centre of the cell
  ! (this is only true if the number of quadrature points is odd)
  qv0 = (nqp_v+1)/2
  qh0 = (nqp_h+1)/2

  ! Step 1: Build integrals of monomials over all cells in advection stencils

  ! Loop over all layers
  layer_loop: do k = 0, nlayers-1

    ! Position vector of centre of this cell
    x0 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = map_wx( df ) + k
      x0(:) = x0(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qh0,qv0)
    end do
    z0 = sqrt(x0(1)**2 + x0(2)**2 + x0(3)**2)*spherical + x0(3)*planar

    ! Compute local order, this is at most vertical_order but reduces near the
    ! top and bottom boundary
    order = min(vertical_order, min(2*(k+1), 2*(nlayers-1 - (k-1))))

    ! Number of monomials to use (all polynomials up to total degree of order)
    nmonomial = (order + 1)
    allocate( int_monomial(nmonomial, nmonomial),  &
              inv_int_monomial(nmonomial, nmonomial), &
              beta(nmonomial), delta(nmonomial), monomial(nmonomial) )

    ! Offset for boundary computations when we move from being upwinded
    ! to downwinded or extrapolated reconstructions
    boundary_offset = 0
    if ( order > 0 ) then
      if ( k == 0 )           boundary_offset =  1
      if ( k == nlayers - 1 ) boundary_offset = -1
    end if

    ! Compute the coefficients of each cell in the stencil for
    ! each edge when this is the upwind cell
    face_loop: do face = 1,2

      int_monomial = 0.0_r_def

      ! Loop over all cells in the stencil
      stencil_loop: do stencil = 1, order+1
        offset = smap(stencil,order) + boundary_offset
        area(stencil) = mdw3(map_w3( 1 ) + k + offset)

        ! Integrate monomials over this cell
        quadrature_loop: do qp = 1, nqp_v
          ! First: Compute physical coordinate of each quadrature point
          xq = 0.0_r_def
          do df = 1, ndf_wx
            ijk = map_wx( df ) + k + offset
            xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qh0,qp)
          end do
          zq = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)*spherical + xq(3)*planar

          ! Second: Compute the local coordinate of each quadrature point from the
          !         physical coordinate
          xx = zq - z0

          ! Third: Compute each needed monomial in terms of the local coordinate
          !        on each quadrature point
          ! Loop over monomials
          do m = 1, nmonomial
            fn = xx**(m-1)
            int_monomial(stencil,m) = int_monomial(stencil,m) &
                                    + wqp_v(qp)*fn*area(stencil)
          end do
        end do quadrature_loop
      end do stencil_loop

      ! Manipulate the integrals of monomials,
      call matrix_invert(int_monomial,inv_int_monomial,nmonomial)

      ! Initialise polynomial coeffficients to zero
      coeff(:,face,map_w3(1)+k) = 0.0_r_def
      ! Loop over quadrature points on this face
      face_quadrature_loop: do qp = 1,nqp_f

        ! Obtain physical coordinates of gauss points on this face
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = map_wx( df ) + k
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*face_basis_wx(1,df,qp,face)
        end do
        zq = sqrt(xq(1)**2 + xq(2)**2 + xq(3)**2)*spherical + xq(3)*planar

        ! Finally obtain local coordinate
        xx = zq - z0

        ! Evaluate polynomial fit
        ! Loop over monomials
        do stencil = 1, order+1
          monomial(stencil) = xx**(stencil-1)
        end do
        do stencil = 1, order+1
          delta(:) = 0.0_r_def
          delta(stencil) = 1.0_r_def
          beta = matmul(inv_int_monomial,delta)
          coeff(stencil,face,map_w3(1)+k) = dot_product(monomial,beta)*area(stencil)
        end do
      end do face_quadrature_loop
    end do face_loop
    deallocate( int_monomial, inv_int_monomial, beta, delta, monomial )
  end do layer_loop

  deallocate( smap )

end subroutine poly1d_vert_flux_coeffs_code

end module poly1d_vert_flux_coeffs_kernel_mod
