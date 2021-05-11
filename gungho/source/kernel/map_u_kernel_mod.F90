!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the rhs for the mapping of the wind field.
!>
!> The kernel computes a very crude approximation to the rhs of the equation u = u0
!> where u0 is the physical wind field. The computational wind field is
!> projected onto using Galerkin projection.
!>
module map_u_kernel_mod

  use argument_mod,            only : arg_type, func_type,       &
                                      GH_FIELD, GH_REAL,         &
                                      GH_INC, GH_READ,           &
                                      ANY_SPACE_9,               &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      GH_BASIS, GH_DIFF_BASIS,   &
                                      CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only : r_def, i_def
  use fs_continuity_mod,       only : W2, W3, Wtheta
  use kernel_mod,              only : kernel_type
  use log_mod,                 only : log_event, LOG_LEVEL_ERROR

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: map_u_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, Wtheta),                   &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),              &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3) &
         /)
    type(func_type) :: meta_funcs(4) = (/                                  &
         func_type(W2,          GH_BASIS),                                 &
         func_type(W3,          GH_BASIS),                                 &
         func_type(Wtheta,      GH_BASIS),                                 &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: map_u_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: map_u_code

contains

!> @brief Compute the right hand side to map the wind field.
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs Right hand side field to compute
!! @param[in] u_lon Longitudinal component of wind in W3
!! @param[in] u_lat Latitudinal component of wind in W3
!! @param[in] u_up Vertical component of wind in Wtheta
!! @param[in] chi_sph_1 1st coordinate in spherical Wchi
!! @param[in] chi_sph_2 2nd coordinate in spherical Wchi
!! @param[in] chi_sph_3 3rd coordinate in spherical Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] ndf_w2 Number of degrees of freedom per cell for w2
!! @param[in] undf_w2 Number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] basis_w2 W2 basis functions evaluated at quadrature points
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] basis_w3 W3 basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_wth Number of degrees of freedom per cell for Wtheta
!! @param[in] undf_wth Number of unique degrees of freedom  for Wtheta
!! @param[in] map_wth Dofmap for the cell at the base of the column for Wtheta
!! @param[in] basis_wt Wtheta basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi_sph Number of degrees of freedom per cell for spherical chi
!! @param[in] undf_chi_sph Number of unique degrees of freedom for spherical chi
!! @param[in] map_chi_sph Dofmap for the cell at the base of the column for spherical chi
!! @param[in] chi_sph_basis Basis functions for spherical Wchi evaluated at
!!                          gaussian quadrature points
!! @param[in] chi_sph_diff_basis Differential of the spherical Wchi basis functions
!!                               evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine map_u_code(nlayers,                                   &
                      rhs,                                       &
                      u_lon, u_lat, u_up,                        &
                      chi_sph_1, chi_sph_2, chi_sph_3, panel_id, &
                      ndf_w2, undf_w2, map_w2, basis_w2,         &
                      ndf_w3, undf_w3, map_w3, basis_w3,         &
                      ndf_wth, undf_wth, map_wth, basis_wt,      &
                      ndf_chi_sph, undf_chi_sph, map_chi_sph,    &
                      chi_sph_basis, chi_sph_diff_basis,         &
                      ndf_pid, undf_pid, map_pid,                &
                      nqp_h, nqp_v, wqp_h, wqp_v                 &
                      )

  use base_mesh_config_mod,       only : geometry,           &
                                         geometry_spherical, &
                                         geometry_planar,    &
                                         topology,           &
                                         topology_fully_periodic
  use chi_transform_mod,          only : chi2llr
  use coordinate_jacobian_mod,    only : coordinate_jacobian
  use coord_transform_mod,        only : sphere2cart_vector

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_chi_sph
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_wth, ndf_w3, ndf_pid
  integer(kind=i_def), intent(in) :: undf_w2, undf_wth, undf_w3, undf_pid
  integer(kind=i_def), intent(in) :: undf_chi_sph
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf_w2),      intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_chi_sph), intent(in) :: map_chi_sph
  integer(kind=i_def), dimension(ndf_w3),      intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wth),     intent(in) :: map_wth
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(3,ndf_w2, nqp_h,nqp_v)     :: basis_w2
  real(kind=r_def), intent(in), dimension(1,ndf_w3, nqp_h,nqp_v)     :: basis_w3
  real(kind=r_def), intent(in), dimension(1,ndf_wth, nqp_h,nqp_v)    :: basis_wt
  real(kind=r_def), intent(in), dimension(1,ndf_chi_sph,nqp_h,nqp_v) :: chi_sph_basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi_sph,nqp_h,nqp_v) :: chi_sph_diff_basis

  real(kind=r_def), dimension(undf_w2),   intent(inout) :: rhs
  real(kind=r_def), dimension(undf_w3),   intent(in)    :: u_lat, u_lon
  real(kind=r_def), dimension(undf_wth),  intent(in)    :: u_up
  real(kind=r_def), dimension(undf_pid),     intent(in) :: panel_id
  real(kind=r_def), dimension(undf_chi_sph), intent(in) :: chi_sph_1, chi_sph_2, chi_sph_3

  real(kind=r_def), dimension(nqp_h), intent(in) ::  wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) ::  wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, k, qp1, qp2, ipanel
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi_sph)     :: chi_sph_1_cell, chi_sph_2_cell, chi_sph_3_cell
  real(kind=r_def), dimension(3)               :: u_physical, u_spherical, coords, llr
  real(kind=r_def)                             :: integrand

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1

    do df = 1, ndf_chi_sph
      chi_sph_1_cell(df) = chi_sph_1( map_chi_sph(df) + k )
      chi_sph_2_cell(df) = chi_sph_2( map_chi_sph(df) + k )
      chi_sph_3_cell(df) = chi_sph_3( map_chi_sph(df) + k )
    end do

    call coordinate_jacobian(ndf_chi_sph,        &
                             nqp_h,              &
                             nqp_v,              &
                             chi_sph_1_cell,     &
                             chi_sph_2_cell,     &
                             chi_sph_3_cell,     &
                             ipanel,             &
                             chi_sph_basis,      &
                             chi_sph_diff_basis, &
                             jacobian,           &
                             dj)
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h

        ! Compute analytical vector wind in physical space
        if ( geometry == geometry_spherical &
             .and. topology == topology_fully_periodic ) then
          ! Need position vector for obtaining (X,Y,Z) components of u_physical
          coords(:) = 0.0_r_def
          do df = 1, ndf_chi_sph
            coords(1) = coords(1) + chi_sph_1_cell(df)*chi_sph_basis(1,df,qp1,qp2)
            coords(2) = coords(2) + chi_sph_2_cell(df)*chi_sph_basis(1,df,qp1,qp2)
            coords(3) = coords(3) + chi_sph_3_cell(df)*chi_sph_basis(1,df,qp1,qp2)
          end do

          call chi2llr(coords(1), coords(2), coords(3), &
                       ipanel, llr(1), llr(2), llr(3))

          u_spherical(1) = u_lon(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_spherical(2) = u_lat(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_spherical(3) = u_up(map_wth(1)+k)*basis_wt(1,1,qp1,qp2) &
                         + u_up(map_wth(2)+k)*basis_wt(1,2,qp1,qp2)
          u_physical     = sphere2cart_vector(u_spherical,llr)

        else if ( geometry == geometry_planar ) then
          u_physical(1) = u_lon(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_physical(2) = u_lat(map_w3(1)+k)*basis_w3(1,1,qp1,qp2)
          u_physical(3) = u_up(map_wth(1)+k)*basis_wt(1,1,qp1,qp2) &
                        + u_up(map_wth(2)+k)*basis_wt(1,2,qp1,qp2)
        else

          call log_event('map_u_kernel is not implemented ' //    &
                         'with your geometry and topology',       &
                         LOG_LEVEL_ERROR)

        end if

        do df = 1, ndf_w2
          integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                         basis_w2(:,df,qp1,qp2)),u_physical)
          rhs(map_w2(df) + k) = rhs(map_w2(df) + k) &
                               + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do

end subroutine map_u_code

end module map_u_kernel_mod
