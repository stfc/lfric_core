!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Linear reconstruction of a scalar field from a coarse to fine mesh
!> @details Performs linear reconstruction of a scalar field (W3 or Wtheta) from
!!          a coarse mesh to a fine mesh, using a stencil to build up the values.
!!          Currently this is hard-wired to use a 2D cross stencil of extent 1,
!!          and when we are at the edge of the domain we just use a constant
!!          reconstruction.
!!          This kernel only works for the lowest-order W3 and Wtheta spaces

module prolong_scalar_linear_kernel_mod

use constants_mod,           only: i_def, r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   GH_COARSE, GH_FINE,        &
                                   STENCIL, CROSS, CELL_COLUMN

implicit none

private

type, public, extends(kernel_type) :: prolong_scalar_linear_kernel_type
   private
   type(arg_type) :: meta_args(2) = (/                                         &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1,       &
                                                          mesh_arg=GH_FINE),   &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2,       &
                                         STENCIL(CROSS),  mesh_arg=GH_COARSE ) &
        /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: prolong_scalar_linear_kernel_code
end type prolong_scalar_linear_kernel_type

public :: prolong_scalar_linear_kernel_code

contains

  !> @brief Linear reconstruction of a scalar field from a coarse to fine mesh
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] fine_field               The fine grid field to write to
  !> @param[in]     coarse_field             Coarse grid field to prolong
  !> @param[in]     stencil_size             Number of cells covered by the stencil
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  subroutine prolong_scalar_linear_kernel_code( nlayers,                 &
                                                cell_map,                &
                                                ncell_fine_per_coarse_x, &
                                                ncell_fine_per_coarse_y, &
                                                ncell_fine,              &
                                                fine_field,              &
                                                coarse_field,            &
                                                stencil_size,            &
                                                stencil_map,             &
                                                ndf,                     &
                                                undf_fine,               &
                                                map_fine,                &
                                                undf_coarse,             &
                                                map_coarse)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: stencil_size
    integer(kind=i_def), intent(in)    :: stencil_map(ndf, stencil_size)
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    real(kind=r_def),    intent(in)    :: coarse_field(undf_coarse)
    real(kind=r_def),    intent(inout) :: fine_field(undf_fine)

    real(kind=r_def)    :: nx, ny, x_fine, y_fine, coeff_a, coeff_b, coeff_c
    integer(kind=i_def) :: df, k, x_idx, y_idx, top_df

    ! Strategy is to approximate the field as:
    ! A*x + B*y + C
    ! between the DoFs of neighbouring coarse cells.
    ! Different coefficients are used for different parts of the coarse cell.
    ! Take the relative positions of DoFs as (xc,yc), (x1,yc) and (xc,y2) for
    ! the points in the coarse cells to use for the interpolation.
    ! The values of the fields at these points are fc, f1 and f2.
    ! Using this we obtain:
    ! A = (fc - f1) / (xc - x1), B = (fc - f2) / (yc - y1)
    ! C = fc - A*xc - B*yc
    ! But now define the coordinate system so that (xc,yc) = (0,0), then we get
    ! A = -(fc - f1) / x1, B = -(fc - f2) / y2, C = fc
    ! These are then used to determine the value at a fine cell centre

    df = 1
    nx = real(ncell_fine_per_coarse_x, r_def)
    ny = real(ncell_fine_per_coarse_y, r_def)

    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf

    ! Only do reconstruction when we have the full stencil
    ! We wouldn't have the full stencil if we were at the edge of a domain,
    ! but for now let's not bother doing any reconstruction there as it is
    ! difficult to work out which edge we're at (and hence which part of the
    ! stencil to use)
    if (stencil_size == 5) then

      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          do k = 0, top_df

            ! C coefficient is the same for all fine cells
            coeff_c = coarse_field(map_coarse(df)+k)

            ! Get relative coordinates of this fine cell within the coarse cell
            x_fine = real(x_idx, r_def) - 0.5_r_def - 0.5_r_def*nx
            y_fine = real(y_idx, r_def) - 0.5_r_def - 0.5_r_def*ny

            ! ---------------------------------------------------------------- !
            ! Find coefficients A and B, depending on location of fine cell
            ! ---------------------------------------------------------------- !
            if (2*x_idx < ncell_fine_per_coarse_x + 1) then
              ! West side of coarse cell
              coeff_a = -(coarse_field(stencil_map(df,2)+k) - coeff_c) / nx
            else if (2*x_idx == ncell_fine_per_coarse_x + 1) then
              ! Central column of coarse cell
              coeff_a = 0.0_r_def
            else
              ! East side of coarse cell
              coeff_a = (coarse_field(stencil_map(df,4)+k) - coeff_c) / nx
            end if

            if (2*y_idx < ncell_fine_per_coarse_y + 1) then
              ! North side of coarse cell
              coeff_b = -(coarse_field(stencil_map(df,5)+k) - coeff_c) / ny
            else if (2*y_idx == ncell_fine_per_coarse_y) then
              ! Central row of coarse cell
              coeff_b = 0.0_r_def
            else
              ! South side of coarse cell
              coeff_b = (coarse_field(stencil_map(df,3)+k) - coeff_c) / ny
            end if

            ! ---------------------------------------------------------------- !
            ! Use coefficients to determine value of fine field
            ! ---------------------------------------------------------------- !
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) = &
              coeff_a*x_fine + coeff_b*y_fine + coeff_c

          end do
        end do
      end do

    else
      ! We're at the edge of the domain. Set all fine values to be the coarse
      ! If in a LAM, these values will be overwritten by LBCs anyway
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          do k = 0, top_df
            ! C coefficient is the same for all fine cells
            coeff_c = coarse_field(map_coarse(df)+k)
            fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) = coeff_c
          end do
        end do
      end do
    end if

  end subroutine prolong_scalar_linear_kernel_code

end module prolong_scalar_linear_kernel_mod
