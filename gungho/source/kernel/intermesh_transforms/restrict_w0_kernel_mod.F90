!!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Restricts a W0 field from a fine mesh to a coarse mesh
!> @details Restrict a W0 field from a fine mesh to coarse mesh. The coarse
!!          field is obtained by taking the values at the coarse cell vertices
!!          from the corresponding vertices on the fine mesh.
!!          Only works for the lowest order elements on quadrilateral cells.
module restrict_w0_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   CELL_COLUMN
use constants_mod,           only: i_def, r_def
use fs_continuity_mod,       only: W0
use kernel_mod,              only: kernel_type
use reference_element_mod,   only: NWB, SWB, SEB, NEB

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: restrict_w0_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                                        &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_FINE  )  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: restrict_w0_kernel_code
end type restrict_w0_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: restrict_w0_kernel_code

contains

  !> @brief Restricts a W0 field from a fine mesh to a coarse mesh
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid W0 field to compute
  !> @param[in]     fine_field               Fine grid W0 field to restrict
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  subroutine restrict_w0_kernel_code(nlayers,                 &
                                     cell_map,                &
                                     ncell_fine_per_coarse_x, &
                                     ncell_fine_per_coarse_y, &
                                     ncell_fine,              &
                                     coarse_field,            &
                                     fine_field,              &
                                     undf_coarse,             &
                                     map_coarse,              &
                                     ndf,                     &
                                     undf_fine,               &
                                     map_fine                 )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse

    ! Fields
    real(kind=r_def),    intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_def),    intent(in)    :: fine_field(undf_fine)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal variables
    integer(kind=i_def) :: k, x_idx, y_idx

    !---------------------------------------------------------------------------
    ! Horizontal components
    !---------------------------------------------------------------------------

    ! Loop from bottom level to top
    do k = 0, nlayers

      ! The rows and columns forming the cell map match the arrangment
      ! of fine cells within the coarse cell

      ! These are aligned as follows with the LFRic directions:
      !         N
      !   |--------------|
      !   |    row 1     |
      !   |c            c|
      !   |o            o|
      ! W |l            l| E
      !   |              |
      !   |1           nx|
      !   |    row ny    |
      !   |--------------|
      !          S

      ! North-west corner. First row, first column of cell map
      x_idx = 1
      y_idx = 1
      coarse_field(map_coarse(NWB)+k) = fine_field(map_fine(NWB,cell_map(x_idx,y_idx))+k)

      ! South-west corner. Last row, first column of cell map
      x_idx = 1
      y_idx = ncell_fine_per_coarse_y
      coarse_field(map_coarse(SWB)+k) = fine_field(map_fine(SWB,cell_map(x_idx,y_idx))+k)

      ! North-east corner. First row, last column of cell map
      x_idx = ncell_fine_per_coarse_x
      y_idx = 1
      coarse_field(map_coarse(NEB)+k) = fine_field(map_fine(NEB,cell_map(x_idx,y_idx))+k)

      ! South-east corner. Last row, last column of cell map
      x_idx = ncell_fine_per_coarse_x
      y_idx = ncell_fine_per_coarse_y
      coarse_field(map_coarse(SEB)+k) = fine_field(map_fine(SEB,cell_map(x_idx,y_idx))+k)

    end do

  end subroutine restrict_w0_kernel_code

end module restrict_w0_kernel_mod
