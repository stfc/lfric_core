!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Calculates a scaling to avoid negative values from intermesh mapping
!> @details Several intergrid mapping processes can generate unwanted negative
!!          values. This can be countered by performing a linear scaling with
!!          a field that is guaranteed not to be negative:
!!          new_field = (1-lambda)*field_to_adjust + lambda*non_negative_field,
!!          where lambda is the scaling factor. Lambda is a field on a coarse
!!          mesh, so that the same correction is applied to all fine cells
!!          within a coarse mesh.
!!          This kernel calculates lambda, as:
!!          lambda = max(-field_to_adjust / (non_negative_field - field_to_adjust))
!!          where the maximum is taken over all fine cells within a coarse cell.
!!          The lambda field is incremented, as it may be desirable to use the
!!          same lambda field for multiple variables (for instance to preserve
!!          linear correlation of mixing ratios).
!!          This kernel is designed for the lowest-order W3 and Wtheta spaces.

module adjust_lambda_kernel_mod

use constants_mod,           only: i_def, r_def
use kernel_mod,              only: kernel_type
use argument_mod,            only: arg_type, GH_REAL,         &
                                   GH_FIELD, CELL_COLUMN,     &
                                   GH_READ, GH_READWRITE,     &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2
use log_mod,                 only: log_event,          &
                                   log_scratch_space,  &
                                   LOG_LEVEL_ALWAYS,   &
                                   LOG_LEVEL_ERROR
implicit none

private

type, public, extends(kernel_type) :: adjust_lambda_kernel_type
   private
   type(arg_type) :: meta_args(3) = (/                                                           &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1, mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_FINE),        &
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_FINE)         &
       /)
  integer :: iterates_over = CELL_COLUMN
contains
  procedure, nopass :: adjust_lambda_kernel_code
end type adjust_lambda_kernel_type

public :: adjust_lambda_kernel_code

contains

  !> @brief Calculates a scaling to avoid negative values from intermesh mapping
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] lambda                   Scaling field on coarse mesh
  !> @param[in]     field_to_adjust          Field on fine mesh, which may have
  !!                                         negative values and will eventually
  !!                                         be corrected.
  !> @param[in]     non_negative_field       Field on fine mesh assured to not
  !!                                         have negative values.
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  subroutine adjust_lambda_kernel_code(nlayers,                 &
                                       cell_map,                &
                                       ncell_fine_per_coarse_x, &
                                       ncell_fine_per_coarse_y, &
                                       ncell_fine,              &
                                       lambda_field,            &
                                       field_to_adjust,         &
                                       non_negative_field,      &
                                       undf_coarse,             &
                                       map_coarse,              &
                                       ndf,                     &
                                       undf_fine,               &
                                       map_fine )

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf
    integer(kind=i_def), intent(in)    :: undf_fine, undf_coarse
    integer(kind=i_def), intent(in)    :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf)
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)
    real(kind=r_def),    intent(inout) :: lambda_field(undf_coarse)
    real(kind=r_def),    intent(in)    :: field_to_adjust(undf_fine)
    real(kind=r_def),    intent(in)    :: non_negative_field(undf_fine)

    real(kind=r_def)    :: pos_value, neg_value, max_lambda, proposed_lambda
    integer(kind=i_def) :: df, top_df, k, x_idx, y_idx

    df = 1

    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_df = nlayers - 2 + ndf

    do k = 0, top_df

      ! Calculate max lambda for the cell to give a positive value
      max_lambda = lambda_field(map_coarse(df)+k)
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          neg_value = field_to_adjust(map_fine(df,cell_map(x_idx,y_idx))+k)
          if (neg_value < 0.0_r_def) then
            pos_value = non_negative_field(map_fine(df,cell_map(x_idx,y_idx))+k)
            proposed_lambda = neg_value / (neg_value - pos_value)
            max_lambda = max(max_lambda, proposed_lambda)
          end if
        end do
      end do

      ! Don't allow values greater than 1
      max_lambda = min(max_lambda, 1.0_r_def)

      ! Only increment lambda if we have exceeded existing value
      lambda_field(map_coarse(df)+k) = max(lambda_field(map_coarse(df)+k), max_lambda)

    end do

  end subroutine adjust_lambda_kernel_code

end module adjust_lambda_kernel_mod
