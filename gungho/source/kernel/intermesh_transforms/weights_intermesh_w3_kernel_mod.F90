!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes weights for conservative prolongation/restriction in W3.
!> @details Computes weights for performing the prolongation operation to take a
!!          W3 field from a coarse mesh to a fine mesh, such that the mass of
!!          the field is conserved. The reciprocal of the weights can be used
!!          for conservative restriction. The weights are stored in a fine W3
!!          field.
!!          This is only designed to work with the lowest order W3 fields.

module weights_intermesh_w3_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_DISCONTINUOUS_SPACE_3, &
                                   CELL_COLUMN
use constants_mod,           only: i_def, r_def
use fs_continuity_mod,       only: W3
use kernel_mod,              only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: weights_intermesh_w3_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3, mesh_arg=GH_FINE ),    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3, mesh_arg=GH_FINE ),    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3, &
                                                  mesh_arg=GH_COARSE )  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: weights_intermesh_w3_kernel_code
end type weights_intermesh_w3_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: weights_intermesh_w3_kernel_code

contains

  !> @brief Computes weights for conservative W3 intermesh mapping.
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] weights_fine             Field containing weights for
  !                                          conservative mapping of W3 fields
  !> @param[in]     mm_w3_fine               Diagonal W3 mass matrix on the fine
  !!                                         mesh (corresponds to cell volumes)
  !> @param[in]     mm_w3_coarse             Diagonal W3 mass matrix on the
  !!                                         coarse mesh
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  subroutine weights_intermesh_w3_kernel_code(nlayers,                 &
                                              cell_map,                &
                                              ncell_fine_per_coarse_x, &
                                              ncell_fine_per_coarse_y, &
                                              ncell_fine,              &
                                              weights_fine,            &
                                              mm_w3_fine,              &
                                              mm_w3_coarse,            &
                                              ndf,                     &
                                              undf_fine,               &
                                              map_fine,                &
                                              undf_coarse,             &
                                              map_coarse               )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=r_def),    intent(inout) :: weights_fine(undf_fine)
    real(kind=r_def),    intent(in)    :: mm_w3_fine(undf_fine)
    real(kind=r_def),    intent(in)    :: mm_w3_coarse(undf_coarse)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: df, k, x_idx, y_idx
    real(kind=r_def)    :: denom

    denom = real(ncell_fine_per_coarse_x*ncell_fine_per_coarse_y, r_def)

    ! Assume lowest-order, so only 1 dof per cell
    df = 1
    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        do k = 0, nlayers-1
          weights_fine(map_fine(df,cell_map(x_idx,y_idx))+k) =                 &
            mm_w3_coarse(map_coarse(df)+k) / denom /                           &
            mm_w3_fine(map_fine(df,cell_map(x_idx,y_idx))+k)
        end do
      end do
    end do

  end subroutine weights_intermesh_w3_kernel_code

end module weights_intermesh_w3_kernel_mod
