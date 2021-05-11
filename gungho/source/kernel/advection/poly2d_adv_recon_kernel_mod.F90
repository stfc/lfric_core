!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes horizontal tracer values through fitting a
!>        high order 2D upwind reconstruction.
!> @details Computes the reconstruction for a tracer field using a high order
!>          polynomial fit to the integrated tracer values. The stencil used
!>          for the polynomial is centred on the upwind cell.
!>          A 2D region stencil centred around the cell is used to build the
!>          polynomial representation of the tracer field.
!>          This method is only valid for lowest order elements.
module poly2d_adv_recon_kernel_mod

use argument_mod,      only : arg_type, func_type,         &
                              reference_element_data_type, &
                              GH_FIELD, GH_SCALAR,         &
                              GH_REAL, GH_INTEGER,         &
                              GH_INC, GH_READ,             &
                              STENCIL, REGION, GH_BASIS,   &
                              CELL_COLUMN, GH_EVALUATOR,   &
                              outward_normals_to_horizontal_faces
use constants_mod,     only : r_def, i_def
use fs_continuity_mod, only : W1, W2, Wtheta
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: poly2d_adv_recon_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                     &
       arg_type(GH_FIELD,  GH_REAL,    GH_INC,  W1),                      &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ, W2),                      &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ, Wtheta, STENCIL(REGION)), &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ, Wtheta),                  &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                          &
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                           &
       /)
  type(func_type) :: meta_funcs(1) = (/                                   &
       func_type(W2, GH_BASIS)                                            &
       /)
  type(reference_element_data_type) :: meta_reference_element(1) = (/     &
       reference_element_data_type( outward_normals_to_horizontal_faces ) &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: poly2d_adv_recon_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: poly2d_adv_recon_code

contains

!> @brief Computes the horizontal polynomial interpolation of a tracer.
!! @param[in]  nlayers Number of layers
!! @param[in,out] reconstruction Reconstructed tracer field to compute
!! @param[in]  wind Wind field
!! @param[in]  tracer Pointwise tracer field to reconstruct
!! @param[in]  coeff Array of polynomial coefficients for interpolation
!! @param[in]  ndf_w1 Number of degrees of freedom per cell
!! @param[in]  undf_w1 Number of unique degrees of freedom for the
!!                     reconstructed field
!! @param[in]  map_w1 Dofmap for the cell at the base of the column
!! @param[in]  ndf_w2 Number of degrees of freedom per cell
!! @param[in]  undf_w2 Number of unique degrees of freedom for the wind field
!! @param[in]  map_w2 Dofmap for the cell at the base of the column
!! @param[in]  basis_w2 Basis function array evaluated at w2 nodes
!! @param[in]  ndf_wt Number of degrees of freedom per cell
!! @param[in]  undf_wt Number of unique degrees of freedom for the tracer field
!! @param[in]  stencil_size Size of the stencil (number of cells)
!! @param[in]  stencil_map Dofmaps for the stencil
!! @param[in]  cells_in_stencil Number of cells needed to compute the polynomial
!!                              fit (may be less than stencil_size)
!! @param[in]  nfaces_re_h Number of horizontal neighbours
!! @param[in]  outward_normals_to_horizontal_faces Vector of normals to the
!!                                                 reference element horizontal
!!                                                 "outward faces"
subroutine poly2d_adv_recon_code( nlayers,              &
                                  reconstruction,       &
                                  wind,                 &
                                  tracer,               &
                                  coeff,                &
                                  ndf_w1,               &
                                  undf_w1,              &
                                  map_w1,               &
                                  ndf_w2,               &
                                  undf_w2,              &
                                  map_w2,               &
                                  basis_w2,             &
                                  ndf_wt,               &
                                  undf_wt,              &
                                  stencil_size,         &
                                  stencil_map,          &
                                  cells_in_stencil,     &
                                  nfaces_re_h,          &
                                  outward_normals_to_horizontal_faces )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                    :: nlayers
  integer(kind=i_def), intent(in)                    :: ndf_wt
  integer(kind=i_def), intent(in)                    :: undf_wt
  integer(kind=i_def), intent(in)                    :: ndf_w1
  integer(kind=i_def), intent(in)                    :: undf_w1
  integer(kind=i_def), dimension(ndf_w1), intent(in) :: map_w1
  integer(kind=i_def), intent(in)                    :: ndf_w2
  integer(kind=i_def), intent(in)                    :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), intent(in)                    :: cells_in_stencil
  integer(kind=i_def), intent(in)                    :: stencil_size
  integer(kind=i_def), intent(in)                    :: nfaces_re_h

  real(kind=r_def), dimension(undf_w1), intent(inout):: reconstruction
  real(kind=r_def), dimension(undf_w2), intent(in)   :: wind
  real(kind=r_def), dimension(undf_wt), intent(in)   :: tracer

  real(kind=r_def), dimension(stencil_size, nfaces_re_h, undf_wt), intent(in) :: coeff

  real(kind=r_def), dimension(3,ndf_w2,ndf_w1), intent(in) :: basis_w2

  integer(kind=i_def), dimension(ndf_wt,stencil_size), intent(in) :: stencil_map

  real(kind=r_def), intent(in)  :: outward_normals_to_horizontal_faces(:,:)

  ! Internal variables
  integer(kind=i_def)                      :: k, df, ij, p
  real(kind=r_def)                         :: direction
  real(kind=r_def), dimension(nfaces_re_h) :: v_dot_n
  real(kind=r_def)                         :: polynomial_tracer

  do df = 1,nfaces_re_h
    v_dot_n(df) = dot_product(basis_w2(:,df,df),outward_normals_to_horizontal_faces(:,df))
  end do

  ij = stencil_map(1,1)

  ! Horizontal tracer reconstruction
  ! Bottom point
  k = 0
  do df = 1,nfaces_re_h
    ! Check if this is the upwind cell
    direction = wind(map_w2(df) + k )*v_dot_n(df)
    if ( direction > 0.0_r_def ) then
      polynomial_tracer = 0.0_r_def
      do p = 1,cells_in_stencil
        polynomial_tracer = polynomial_tracer &
                          + tracer( stencil_map(1,p) + k )*coeff( p, df, ij+k )
      end do
      reconstruction(map_w1(df) + k ) = polynomial_tracer
    end if
  end do

  do k = 1, nlayers - 1
    do df = 1,nfaces_re_h
      ! Check if this is the upwind cell
      direction = (wind(map_w2(df) + k ) + wind(map_w2(df) + k-1 ))*v_dot_n(df)
      if ( direction > 0.0_r_def ) then
        polynomial_tracer = 0.0_r_def
        do p = 1,cells_in_stencil
          polynomial_tracer = polynomial_tracer &
                            + tracer( stencil_map(1,p) + k )*coeff( p, df, ij+k )
        end do
        reconstruction(map_w1(df) + k ) = polynomial_tracer
      end if
    end do
  end do
  ! Final point
  k = nlayers - 1
  do df = 1,nfaces_re_h
    ! Check if this is the upwind cell
    direction = wind(map_w2(df) + k )*v_dot_n(df)
    if ( direction > 0.0_r_def ) then
      polynomial_tracer = 0.0_r_def
      do p = 1,cells_in_stencil
        polynomial_tracer = polynomial_tracer &
                          + tracer( stencil_map(2,p) + k )*coeff( p, df, ij+k+1 )
      end do
      reconstruction(map_w1(df) + k + 1 ) = polynomial_tracer
    end if
  end do

end subroutine poly2d_adv_recon_code

end module poly2d_adv_recon_kernel_mod
