!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Compute the coefficients for reconstructing a
!>        a  2D horizontal upwind polynomial representation of a tracer field on the
!>        faces of a cell
!> @details Compute the coefficients of the flux of a tracer density field using a high order
!>          2D polynomial fit to the integrated tracer values over a given stencil.
!>          The stencil used for the polynomial is centred on the upwind cell for each edge
!>          A symmetric polynomial is used containing all monomials up to the
!>          desired order, i.e. order = 2: 1 + x + y + x^2 + xy + y^2
!>          This is exactly fitted over the central cell and then in a least
!>          squares manner over the remaining cells in the stencil
!>          The methodology closely follows that of Thuburn et.al GMD 2014
!>          This method is only valid for lowest order elements
module poly2d_flux_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type,                  &
                              GH_FIELD, GH_INTEGER, GH_REAL,        &
                              GH_WRITE, GH_READ,                    &
                              ANY_SPACE_1,                          &
                              ANY_DISCONTINUOUS_SPACE_3,            &
                              GH_BASIS, CELLS, GH_QUADRATURE_XYoZ,  &
                              GH_QUADRATURE_face,                   &
                              STENCIL, REGION
use constants_mod,     only : r_def, i_def, l_def
use fs_continuity_mod, only : W3
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: poly2d_flux_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                                  &
       arg_type(GH_FIELD,   GH_WRITE, W3),                             &
       arg_type(GH_FIELD,   GH_READ,  W3,          STENCIL(REGION)),   &
       arg_type(GH_FIELD*3, GH_READ,  ANY_SPACE_1, STENCIL(REGION)),   &
       arg_type(GH_FIELD,   GH_READ,  ANY_DISCONTINUOUS_SPACE_3, STENCIL(REGION)),   &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_INTEGER, GH_READ),                                  &
       arg_type(GH_REAL,    GH_READ)                                   &
       /)
  type(func_type) :: meta_funcs(1) = (/                                &
       func_type(ANY_SPACE_1, GH_BASIS)                                &
       /)
  integer :: iterates_over = CELLS
  integer :: gh_shape(2) = (/ GH_QUADRATURE_XYoZ, GH_QUADRATURE_face /)
contains
  procedure, nopass :: poly2d_flux_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public poly2d_flux_coeffs_code
contains

!>@brief Compute the coefficients needed for a 2D horizontal reconstruction
!>       of a tracer field on horizontal faces
!>@param[in] nlayers Number of vertical layers
!>@param[out] coeff Array of fields to store the coefficients for the polynomial
!!                  reconstruction
!>@param[in]  mdw3 Mass matrix diagonal for the W3 space, this is used to give
!!                 the cell volume
!>@param[in] chi1 1st component of the physical coordinate field
!>@param[in] chi2 2nd component of the physical coordinate field
!>@param[in] chi3 3rd component of the physical coordinate field
!>@param[in] panel_id Id of the cubed sphere panel for each column
!>@param[in] ndf_w3 Number of degrees of freedom per cell for W3
!>@param[in] undf_w3 Total number of degrees of freedom for W3
!>@param[in] stencil_size_w3 Number of cells in the W3 stencil
!>@param[in] smap_w3 Stencil dofmap of the W3 stencil
!>@param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!>@param[in] undf_wx Total number of degrees of freedom for the coordinate space
!>@param[in] stencil_size_wx Number of cells in the coordinate space stencil
!>@param[in] smap_wx Stencil dofmap of the coordinate space stencil
!>@param[in] basis_wx Basis function of the coordinate space evaluated on
!!                    quadrature points
!>@param[in] face_basis_wx Basis function of the coordinate space evaluated on
!!                         quadrature points on the horizontal faces
!>@param[in] ndf_pid Number of degrees of freedom per cell for the panel_id space
!>@param[in] undf_pid Total number of degrees of freedom per cell for the panel_id space
!>@param[in] stencil_size_pid Number of cells in the stencil for the panel_id space
!>@param[in] smap_pid Stencil dofmap for the panel_id space
!>@param[in] order Polynomial order for flux computations
!>@param[in] nfaces_h Number of horizontal neighbours
!>@param[in] cells_in_stencil Number of cells in the stencil to use for the
!!                            reconstruction in this column (may be smaller than
!!                            stencil_size_w3)
!>@param[in] chi3_min Global minimum of chi3 array
!>@param[in] nqp_h Number of horizontal quadrature points
!>@param[in] nqp_v Number of vertical quadrature points
!>@param[in] wqp_h Weights of horizontal quadrature points
!>@param[in] wqp_v Weights of vertical quadrature points
!>@param[in] nfaces_qr Number of faces in the quadrature rule
!>@param[in] nqp_f Number of face quadrature points
!>@param[in] wqp_f Weights of face quadrature points
subroutine poly2d_flux_coeffs_code(nlayers,                    &
                                   coeff,                      &
                                   mdw3,                       &
                                   chi1, chi2, chi3,           &
                                   panel_id,                   &
                                   ndf_w3,                     &
                                   undf_w3,                    &
                                   stencil_size_w3,            &
                                   smap_w3,                    &
                                   ndf_wx,                     &
                                   undf_wx,                    &
                                   stencil_size_wx,            &
                                   smap_wx,                    &
                                   basis_wx,                   &
                                   face_basis_wx,              &
                                   ndf_pid, undf_pid,          &
                                   stencil_size_pid,           &
                                   smap_pid,                   &
                                   order,                      &
                                   nfaces_h,                   &
                                   cells_in_stencil,           &
                                   chi3_min,                   &
                                   nqp_h, nqp_v, wqp_h, wqp_v, &
                                   nfaces_qr, nqp_f, wqp_f )


  use matrix_invert_mod,         only: matrix_invert
  use cross_product_mod,         only: cross_product
  use base_mesh_config_mod,      only: geometry, &
                                       geometry_spherical
  use poly_helper_functions_mod, only: buildadvcoeff, &
                                       local_distance_2d
  use finite_element_config_mod, only: spherical_coord_system,     &
                                       spherical_coord_system_abh
  use coord_transform_mod,       only: alphabetar2xyz
  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: order, cells_in_stencil
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3, &
                                     ndf_wx, undf_wx, &
                                     ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: nqp_v, nqp_h, nqp_f, nfaces_qr, nfaces_h
  integer(kind=i_def), intent(in) :: stencil_size_w3, stencil_size_wx, stencil_size_pid

  integer(kind=i_def), dimension(ndf_w3, stencil_size_w3),  intent(in) :: smap_w3
  integer(kind=i_def), dimension(ndf_wx, stencil_size_wx),  intent(in) :: smap_wx
  integer(kind=i_def), dimension(ndf_pid,stencil_size_pid), intent(in) :: smap_pid

  real(kind=r_def), dimension(undf_w3),                            intent(in)  :: mdw3
  real(kind=r_def), dimension(undf_wx),                            intent(in)  :: chi1, chi2, chi3
  real(kind=r_def), dimension(stencil_size_w3, nfaces_h, undf_w3), intent(out) :: coeff
  real(kind=r_def), dimension(undf_pid),                           intent(in)  :: panel_id

  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v),     intent(in) :: basis_wx
  real(kind=r_def), dimension(1,ndf_wx,nqp_f,nfaces_qr), intent(in) :: face_basis_wx

  real(kind=r_def), dimension(nqp_h),           intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v),           intent(in) ::  wqp_v
  real(kind=r_def), dimension(nqp_f,nfaces_qr), intent(in) ::  wqp_f

  real(kind=r_def), intent(in) :: chi3_min

  ! Local variables
  logical(kind=l_def) :: spherical
  integer(kind=i_def) :: ispherical, ipanel
  integer(kind=i_def) :: k, ijk, df, qv0, qh0, stencil, nmonomial, qp, &
                         px, py, m, face

  real(kind=r_def)                              :: fn, poly
  real(kind=r_def),              dimension(2)   :: xx
  real(kind=r_def),              dimension(3)   :: x0, x1, xq, xn1, abr, r0
  real(kind=r_def), allocatable, dimension(:,:) :: int_monomial
  real(kind=r_def),              dimension(stencil_size_w3) :: area

  ! Ensure that r0 + r is positive
  r0 = (/ 0.0_r_def, 0.0_r_def, 1.0_r_def + abs(chi3_min) /)

  if ( geometry == geometry_spherical ) then
    spherical = .true.
    ispherical = 1_i_def
  else
    spherical = .false.
    ispherical = 0_i_def
  end if

  area(:) = 1.0_r_def

  ! Avoid compile warning for unused variables
  fn = wqp_v(1)

  ! Number of monomials to use (all polynomials up to total degree of order)
  nmonomial = (order + 1)*(order + 2)/2

  ! Index of quadrature point in the centre of the cell
  ! (this is only true if the number of quadrature points is odd
  qv0 = (nqp_v+1)/2
  qh0 = (nqp_h+1)/2

  ! Step 1: Build integrals of monomials over all cells in advection stencils
  ! Initialize to zero
  allocate( int_monomial(stencil_size_w3, nmonomial) )

  ! Loop over all layers
  layer_loop: do k = 0, nlayers-1

    int_monomial = 0.0_r_def

    ! Position vector of centre of this cell
    x0 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = smap_wx( df, 1) + k
      x0(:) = x0(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qh0,qv0)
    end do
    ! Find direction of first neighbour to establish axes of
    ! Local coordinate system
    ! Position vector of neighbour cell centre
    x1 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = smap_wx( df, 2) + k
      x1(:) = x1(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qh0,qv0)
    end do

    if ( spherical_coord_system == spherical_coord_system_abh ) then
      ! Convert x0 & x1 to XYZ coordinate system
      ipanel = int(panel_id(smap_pid(1,1)), i_def)
      abr = x0 + r0
      call alphabetar2xyz(abr(1), abr(2), abr(3), &
                          ipanel, x0(1), x0(2), x0(3))
      ipanel = int(panel_id(smap_pid(1,2)), i_def)
      abr = x1 + r0
      call alphabetar2xyz(abr(1), abr(2), abr(3), &
                          ipanel, x1(1), x1(2), x1(3))
    end if

    x1(3) = ispherical*x1(3) + (1_i_def - ispherical)*x0(3)
    ! Unit normal to plane containing points 0 and 1
    xn1 = cross_product(x0,x1)
    xn1 = xn1/sqrt(xn1(1)**2 + xn1(2)**2 + xn1(3)**2)

    ! Loop over all cells in stencils
    stencil_loop: do stencil = 1, cells_in_stencil

      ! Integrate monomials over this cell
      area(stencil) = mdw3(smap_w3( 1, stencil) + k)
      quadrature_loop: do qp = 1, nqp_h
        ! First: Compute physical coordinate of each quadrature point
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = smap_wx( df, stencil) + k
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qp,qv0)
        end do
        if ( spherical_coord_system == spherical_coord_system_abh ) then
          ! Convert xq to XYZ coordinate system
          ipanel = int(panel_id(smap_pid(1,stencil)), i_def)
          abr = xq + r0
          call alphabetar2xyz(abr(1), abr(2), abr(3), &
                              ipanel, xq(1), xq(2), xq(3))
        end if

        ! Second: Compute the local coordinate of each quadrature point from the
        !         physical coordinate
        xx = local_distance_2d(x0, xq, xn1, spherical)

        ! Third: Compute each needed monomial in terms of the local coordinate
        !        on each quadrature point
        ! Loop over monomials
        px = 0
        py = 0
        do m = 1, nmonomial
          fn = (xx(1)**px)*(xx(2)**py)
          int_monomial(stencil,m) = int_monomial(stencil,m) &
                                  + wqp_h(qp)*fn*area(stencil)
          px = px - 1
          py = py + 1
          if (px < 0) then
            px = py
            py = 0
          end if
        end do
      end do quadrature_loop
    end do stencil_loop
    ! Manipulate the integrals of monomials
    call buildadvcoeff(int_monomial, stencil_size_w3, nmonomial)

    ! Now compute the coefficients of each cell in the stencil for
    ! each edge when this is the upwind cell
    face_loop: do face = 1,nfaces_h
      ! Initialise polynomial coeffficients to zero
      coeff(:,face,smap_w3(1,1)+k) = 0.0_r_def
      ! Loop over quadrature points on this face
      face_quadrature_loop: do qp = 1,nqp_f

        ! Obtain physical coordinates of gauss points on this face
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = smap_wx( df, 1) + k
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*face_basis_wx(1,df,qp,face)
        end do
        if ( spherical_coord_system == spherical_coord_system_abh ) then
          ! Convert xq to XYZ coordinate system
          ipanel = int(panel_id(smap_pid(1,1)), i_def)
          abr = xq + r0
          call alphabetar2xyz(abr(1), abr(2), abr(3), &
                              ipanel, xq(1), xq(2), xq(3))
        end if

        xx = local_distance_2d(x0, xq, xn1, spherical)

        ! Evaluate polynomial fit
        ! Loop over monomials
        do stencil = 1, cells_in_stencil
          poly = 0.0_r_def
          px = 0
          py = 0
          do m = 1, nmonomial
            fn = (xx(1)**px)*(xx(2)**py)
            poly = poly + int_monomial(stencil,m)*fn
            px = px - 1
            py = py + 1
            if (px < 0) then
              px = py
              py = 0
            end if
          end do
          coeff(stencil,face,smap_w3(1,1)+k) = &
            coeff(stencil,face,smap_w3(1,1)+k) + wqp_f(qp,face)*poly*area(stencil)
        end do
      end do face_quadrature_loop

    end do face_loop
  end do layer_loop

  deallocate( int_monomial )
end subroutine poly2d_flux_coeffs_code

end module poly2d_flux_coeffs_kernel_mod
