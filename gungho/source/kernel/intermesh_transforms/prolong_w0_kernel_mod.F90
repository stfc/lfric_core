!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Prolongs a W0 field from a coarse mesh to a fine mesh
!> @details Prolong a W0 field on a coarse mesh into a W0 field on a fine mesh.
!!          The new values of the W0 fine field which coincide with the vertices
!!          of the coarse cell will have the corresponding values of the coarse
!!          cell. Values interior to the coarse cell are obtained by
!!          interpolating.
!!          Only works for the lowest order elements on quadrilateral cells.
module prolong_w0_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   CELL_COLUMN
use constants_mod,           only: i_def, r_def
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

type, public, extends(kernel_type) :: prolong_w0_kernel_type
  private
  type(arg_type) :: meta_args(2) = (/                                                        &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1, mesh_arg=GH_FINE ),  &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_COARSE ) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: prolong_w0_kernel_code
end type prolong_w0_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: prolong_w0_kernel_code

contains

  !> @brief Prolongs a W0 field from a coarse mesh to a fine mesh
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] fine_field               Fine grid W0 field to restrict
  !> @param[in]     coarse_field             Coarse grid W0 field to compute
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  subroutine prolong_w0_kernel_code(nlayers,                 &
                                    cell_map,                &
                                    ncell_fine_per_coarse_x, &
                                    ncell_fine_per_coarse_y, &
                                    ncell_fine,              &
                                    fine_field,              &
                                    coarse_field,            &
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
    real(kind=r_def),    intent(inout) :: fine_field(undf_fine)
    real(kind=r_def),    intent(in)    :: coarse_field(undf_coarse)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: k, x_idx, y_idx
    real(kind=r_def)    :: weight_NE, weight_NW, weight_SE, weight_SW

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

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x

        ! First job is to orient fine cell within coarse cell
        ! Do SW DoF for every fine cell
        weight_SW = real(y_idx, r_def)                                       &
                    * real(ncell_fine_per_coarse_x - x_idx + 1_i_def, r_def) &
                    / real(ncell_fine_per_coarse_x * ncell_fine_per_coarse_y, r_def)
        weight_NW = real(ncell_fine_per_coarse_y - y_idx, r_def)             &
                    * real(ncell_fine_per_coarse_x - x_idx + 1_i_def, r_def) &
                    / real(ncell_fine_per_coarse_x * ncell_fine_per_coarse_y, r_def)
        weight_NE = real(ncell_fine_per_coarse_y - y_idx, r_def)             &
                    * real(x_idx - 1_i_def, r_def)                           &
                    / real(ncell_fine_per_coarse_x * ncell_fine_per_coarse_y, r_def)
        weight_SE = real(y_idx, r_def) * real(x_idx - 1_i_def, r_def)        &
                    / real(ncell_fine_per_coarse_x * ncell_fine_per_coarse_y, r_def)

        ! Loop from bottom of bottom layer to top of top layer
        do k = 0, nlayers

          fine_field(map_fine(SWB,cell_map(x_idx,y_idx))+k) =             &
                          weight_NE * coarse_field(map_coarse(NEB)+k)     &
                          + weight_SE * coarse_field(map_coarse(SEB)+k)   &
                          + weight_SW * coarse_field(map_coarse(SWB)+k)   &
                          + weight_NW * coarse_field(map_coarse(NWB)+k)
        end do

        ! Do SE DoF when in last column
        if ( x_idx == ncell_fine_per_coarse_x ) then
          ! Weights are just based on SE and NE coarse values
          weight_SE = real(y_idx, r_def)  &
                      / real(ncell_fine_per_coarse_y, r_def)
          weight_NE = 1.0_r_def - weight_SE

          do k = 0, nlayers
            fine_field(map_fine(SEB,cell_map(x_idx,y_idx))+k) =           &
                            weight_NE * coarse_field(map_coarse(NEB)+k)   &
                            + weight_SE * coarse_field(map_coarse(SEB)+k)
          end do
        end if

        ! Do NW DoF when in first row
        if ( y_idx == 1_i_def ) then
          ! Weights are just based on NW and NE coarse values
          weight_NW = real(ncell_fine_per_coarse_x - x_idx + 1_i_def, r_def) &
                      / real(ncell_fine_per_coarse_x, r_def)
          weight_NE = 1.0_r_def - weight_NW

          do k = 0, nlayers
            fine_field(map_fine(NWB,cell_map(x_idx,y_idx))+k) =           &
                            weight_NE * coarse_field(map_coarse(NEB)+k)   &
                            + weight_NW * coarse_field(map_coarse(NWB)+k)
          end do

          ! Do NE DoF when in first row and last column
          if ( x_idx == ncell_fine_per_coarse_x ) then
            do k = 0, nlayers
              ! Value is just NE value of coarse cell
              fine_field(map_fine(NEB,cell_map(x_idx,y_idx))+k) =         &
                            coarse_field(map_coarse(NEB)+k)
            end do
          end if
        end if

      end do ! Loop through rows of cell map
    end do ! Loop through columns of cell map

  end subroutine prolong_w0_kernel_code

end module prolong_w0_kernel_mod
