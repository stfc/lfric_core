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
!> @details Compute the coefficients of the advective update of a tracer field using a high order
!>          2D polynomial fit to the integrated tracer values over a given stencil.
!>          The stencil used for the polynomial is centred on the upwind cell for each edge
!>          A symmetric polynomial is used containing all monomials up to the
!>          desired order, i.e. order = 2: 1 + x + y + x^2 + xy + y^2
!>          This is exactly fitted over all cells in the stencil
!>          The methodology closely follows that of Thuburn et.al GMD 2014 for
!>          2D reconstructions
!>          This method is only valid for lowest order elements
module poly2d_advective_coeffs_kernel_mod

use argument_mod,      only : arg_type, func_type,       &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_WRITE, GH_READ,         &
                              ANY_SPACE_1, CELL_COLUMN,  &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              STENCIL, REGION, GH_BASIS, &
                              GH_QUADRATURE_XYoZ,        &
                              GH_QUADRATURE_edge
use constants_mod,     only : r_def, i_def, l_def
use fs_continuity_mod, only : Wtheta
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly2d_advective_coeffs_kernel_type
  private
  type(arg_type) :: meta_args(8) = (/                                            &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, Wtheta),                       &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  Wtheta,      STENCIL(REGION)), &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  ANY_SPACE_1, STENCIL(REGION)), &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3,     &
                                                  STENCIL(REGION)),              &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                                &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                                &
       arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),                                &
       arg_type(GH_SCALAR,  GH_REAL,    GH_READ)                                 &
       /)
  type(func_type) :: meta_funcs(1) = (/                                          &
       func_type(ANY_SPACE_1, GH_BASIS)                                          &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape(2) = (/ GH_QUADRATURE_XYoZ, GH_QUADRATURE_edge /)
contains
  procedure, nopass :: poly2d_advective_coeffs_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly2d_advective_coeffs_code
contains

!> @brief Compute the coefficients needed for a 1D horizontal reconstruction
!>        of a tracer field on horizontal faces
!> @param[in] nlayers Number of vertical layers
!> @param[in,out] coeff Array of fields to store the coefficients for the
!!                      polynomial reconstruction
!> @param[in] mdwt Mass matrix diagonal for the Wtheta space, this is used to
!!                 give the cell volume
!> @param[in] chi1 1st component of the physical coordinate field
!> @param[in] chi2 2nd component of the physical coordinate field
!> @param[in] chi3 3rd component of the physical coordinate field
!> @param[in] panel_id Id of the cubed sphere panel for each column
!> @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wt Total number of degrees of freedom for Wtheta
!> @param[in] stencil_size_wt Number of cells in the Wtheta stencil
!> @param[in] smap_wt Stencil dofmap of the Wtheta stencil
!> @param[in] ndf_wx Number of degrees of freedom per cell for the coordinate space
!> @param[in] undf_wx Total number of degrees of freedom for the coordinate space
!> @param[in] stencil_size_wx Number of cells in the coordinate space stencil
!> @param[in] smap_wx Stencil dofmap of the coordinate space stencil
!> @param[in] basis_wx Basis function of the coordinate space evaluated on
!!                     quadrature points. The vertical aspect must be on GLL points
!> @param[in] edge_basis_wx Basis function of the coordinate space evaluated on
!!                          quadrature points on the horizontal edges,
!> @param[in] ndf_pid Number of degrees of freedom per cell for the panel_id space
!> @param[in] undf_pid Total number of degrees of freedom per cell for the panel_id space
!> @param[in] stencil_size_pid Number of cells in the stencil for the panel_id space
!> @param[in] smap_pid Stencil dofmap for the panel_id space
!> @param[in] order Polynomial order for flux computations
!> @param[in] nfaces_h Number of horizontal neighbours
!> @param[in] cells_in_stencil Number of cells in the stencil to use for the
!!                             reconstruction in this column (may be smaller than
!!                             stencil_size_wt)
!> @param[in] transform_radius A radius used for transforming to spherically-based
!!                             coords. For Cartesian coordinates this is zero, but
!!                             for spherical coordinates it is the global minimum
!!                             of the height field plus 1.
!> @param[in] nqp_h Number of horizontal quadrature points
!> @param[in] nqp_v Number of vertical quadrature points
!> @param[in] wqp_h Weights of horizontal quadrature points
!> @param[in] wqp_v Weights of vertical quadrature points
!> @param[in] nedges_qr Number of edges in the quadrature rule
!> @param[in] nqp_e Number of edge quadrature points
!> @param[in] wqp_e Weights of edge quadrature points
subroutine poly2d_advective_coeffs_code(nlayers,                    &
                                        coeff,                      &
                                        mdwt,                       &
                                        chi1, chi2, chi3,           &
                                        panel_id,                   &
                                        ndf_wt,                     &
                                        undf_wt,                    &
                                        stencil_size_wt,            &
                                        smap_wt,                    &
                                        ndf_wx,                     &
                                        undf_wx,                    &
                                        stencil_size_wx,            &
                                        smap_wx,                    &
                                        basis_wx,                   &
                                        edge_basis_wx,              &
                                        ndf_pid,                    &
                                        undf_pid,                   &
                                        stencil_size_pid,           &
                                        smap_pid,                   &
                                        order,                      &
                                        nfaces_h,                   &
                                        cells_in_stencil,           &
                                        transform_radius,           &
                                        nqp_h, nqp_v, wqp_h, wqp_v, &
                                        nedges_qr, nqp_e, wqp_e )

  use matrix_invert_mod,         only: matrix_invert
  use cross_product_mod,         only: cross_product
  use base_mesh_config_mod,      only: geometry, &
                                       geometry_spherical
  use poly_helper_functions_mod, only: buildadvcoeff, &
                                       local_distance_2d
  use chi_transform_mod,         only: chir2xyz

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: order, cells_in_stencil, nfaces_h
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt, &
                                     ndf_wx, undf_wx, &
                                     ndf_pid, undf_pid
  integer(kind=i_def), intent(in) :: nqp_v, nqp_h, nqp_e, nedges_qr
  integer(kind=i_def), intent(in) :: stencil_size_wt, stencil_size_wx, stencil_size_pid

  integer(kind=i_def), dimension(ndf_wt, stencil_size_wt),  intent(in) :: smap_wt
  integer(kind=i_def), dimension(ndf_wx, stencil_size_wx),  intent(in) :: smap_wx
  integer(kind=i_def), dimension(ndf_pid,stencil_size_pid), intent(in) :: smap_pid

  real(kind=r_def), dimension(undf_wt),                            intent(in)    :: mdwt
  real(kind=r_def), dimension(undf_wx),                            intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(stencil_size_wt, nfaces_h, undf_wt), intent(inout) :: coeff
  real(kind=r_def), dimension(undf_pid),                           intent(in)    :: panel_id

  real(kind=r_def), dimension(1,ndf_wx,nqp_h,nqp_v),     intent(in) :: basis_wx
  real(kind=r_def), dimension(1,ndf_wx,nqp_e,nedges_qr), intent(in) :: edge_basis_wx

  real(kind=r_def), dimension(nqp_h),           intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v),           intent(in) :: wqp_v
  real(kind=r_def), dimension(nqp_e,nedges_qr), intent(in) :: wqp_e

  real(kind=r_def), intent(in) :: transform_radius

  ! Local variables
  logical(kind=l_def) :: spherical
  integer(kind=i_def) :: ispherical, ipanel
  integer(kind=i_def) :: k, ijk, df, qf0, stencil, nmonomial, qp, &
                         m, edge, px, py
  integer(kind=i_def), dimension(0:nlayers) :: kx, vert_face, qv
  real(kind=r_def)                                          :: fn, poly, z0
  real(kind=r_def),              dimension(2)               :: xx
  real(kind=r_def),              dimension(3)               :: x0, x1, xq, xn1, chi_sph, r0
  real(kind=r_def), allocatable, dimension(:,:)             :: int_monomial
  real(kind=r_def),              dimension(stencil_size_wt) :: area

  ! Radius correction for transforming to (X,Y,Z) coordinates
  r0 = (/ 0.0_r_def, 0.0_r_def, transform_radius /)

  if ( geometry == geometry_spherical ) then
    spherical = .true.
    ispherical = 1_i_def
    z0 = 0.0_r_def
  else
    spherical = .false.
    ispherical = 0_i_def
    z0 = 1.0_r_def
  end if

  area(:) = 1.0_r_def

  ! Avoid compile warning for unused variables
  fn = wqp_v(1)

  ! Number of monomials to use (all polynomials up to total degree of order)
  nmonomial = (order + 1)*(order + 2)/2

  ! Index of quadrature point in the centre of the vertical face
  ! (this is only true if the number of quadrature points is odd
  qf0 = (nqp_h+1)/2

  ! Step 1: Build integrals of monomials over all cells in advection stencils
  ! Initialize to zero
  allocate( int_monomial(stencil_size_wt, nmonomial) )

  ! Set up indexing arrays to account for special nature of
  ! top point
  vert_face = 0
  vert_face(nlayers) = 4
  kx = (/(k, k=0, nlayers)/)
  kx(nlayers) = nlayers-1
  qv = 1
  qv(nlayers) = nqp_v

  ! Loop over all layers: goes to nalyers to pick up top point
  layer_loop: do k = 0, nlayers

    int_monomial = 0.0_r_def

    ! Position vector of bottom of this cell unless very last point in
    ! which case use the top of the cell
    x0 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = smap_wx(df, 1) + kx(k)
      x0(:) = x0(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qf0,qv(k))
    end do

    ! Convert x0 to XYZ coordinate system
    ipanel = int(panel_id(smap_pid(1,1)), i_def)
    chi_sph = x0 + r0
    call chir2xyz(chi_sph(1), chi_sph(2), chi_sph(3), &
                  ipanel, x0(1), x0(2), x0(3))

    ! Avoid issues when x0(3) == 0
    if ( k == 0) x0(3) = x0(3) + z0

    ! Find direction of first neighbour to establish axes of
    ! Local coordinate system
    ! Position vector of neighbour cell centre
    x1 = 0.0_r_def
    do df = 1, ndf_wx
      ijk = smap_wx( df, 2) + kx(k)
      x1(:) = x1(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qf0,qv(k))
    end do

    ! Convert x1 to XYZ coordinate system
    ipanel = int(panel_id(smap_pid(1,2)), i_def)
    chi_sph = x1 + r0
    call chir2xyz(chi_sph(1), chi_sph(2), chi_sph(3), &
                  ipanel, x1(1), x1(2), x1(3))

    x1(3) = ispherical*x1(3) + (1_i_def - ispherical)*x0(3)
    ! Unit normal to plane containing points 0 and 1
    xn1 = cross_product(x0,x1)
    xn1 = xn1/sqrt(xn1(1)**2 + xn1(2)**2 + xn1(3)**2)

    ! Loop over all cells in the stencil
    stencil_loop: do stencil = 1, cells_in_stencil
      area(stencil) = mdwt(smap_wt( 1, stencil) + k)
      if ( k > 0 .and. k < nlayers ) &
        area(stencil) =  area(stencil) &
          + mdwt(smap_wt( 2, stencil) + k-1)
      ! Integrate monomials over this cell
      quadrature_loop: do qp = 1, nqp_h
        ! First: Compute physical coordinate of each quadrature point
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = smap_wx(df, stencil) + kx(k)
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*basis_wx(1,df,qp,qv(k))
        end do

        ! Convert xq to XYZ coordinate system
        ipanel = int(panel_id(smap_pid(1, stencil)), i_def)
        chi_sph = xq + r0
        call chir2xyz(chi_sph(1), chi_sph(2), chi_sph(3), &
                      ipanel, xq(1), xq(2), xq(3))

        ! Avoid issues when x0(3) == 0
        if ( k == 0) xq(3) = xq(3) + z0

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
    call buildadvcoeff(int_monomial, stencil_size_wt, nmonomial)

    ! Now compute the coefficients of each cell in the stencil for
    ! each edge when this is the upwind cell
    edge_loop: do edge = 1,nfaces_h
      ! Initialise polynomial coeffficients to zero
      coeff(:,edge,smap_wt(1,1)+k) = 0.0_r_def
      ! Loop over quadrature points on this edge
      edge_quadrature_loop: do qp = 1,nqp_e

        ! Obtain physical coordinates of gauss points on this edge
        xq = 0.0_r_def
        do df = 1, ndf_wx
          ijk = smap_wx(df, 1) + kx(k)
          xq(:) = xq(:) + (/ chi1(ijk), chi2(ijk), chi3(ijk) /)*edge_basis_wx(1,df,qp,edge + vert_face(k))
        end do

        ! Convert xq to XYZ coordinate system
        ipanel = int(panel_id(smap_pid(1, 1)), i_def)
        chi_sph = xq + r0
        call chir2xyz(chi_sph(1), chi_sph(2), chi_sph(3), &
                      ipanel, xq(1), xq(2), xq(3))

        ! Obtain local coordinates of gauss points on this edge
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
          coeff(stencil,edge,smap_wt(1,1)+k) = &
            coeff(stencil,edge,smap_wt(1,1)+k) + wqp_e(qp,edge)*poly*area(stencil)
        end do
      end do edge_quadrature_loop
    end do edge_loop
  end do layer_loop
  deallocate( int_monomial )

end subroutine poly2d_advective_coeffs_code

end module poly2d_advective_coeffs_kernel_mod
