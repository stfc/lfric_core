!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Perform the restriction operation from a fine grid W2 field to a
!!        coarse grid W2 field, using a mask e.g. a limited area mask, for
!!        an intensive field.
!> @details Restrict the W2 fine grid field over a number of cells into a
!!          W2 coarse grid field. The fine grid cells must be exactly
!!          nested in a coarse grid cell. The coarse field is obtained by
!!          summing contributions from the fine field multiplied by weights.
!!          This method is only designed for the lowest order W2 spaces.
!!          Deal with the N,S,E,W dofs separately. For example, for the W dofs
!!          For a coarse mesh cell j and W dof, that contains Nf fine mesh cells
!!          coarse(j,W) = sum_{i=1,Nf} mask(i,W) * fine(i,W) /  sum_{i=1,Nf} mask(i,W)
module restrict_masked_w2_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_INC,           &
                                   GH_COARSE, GH_FINE,        &
                                   ANY_SPACE_2, CELL_COLUMN
use constants_mod,           only: i_def, r_def
use fs_continuity_mod,       only: W2
use kernel_mod,              only: kernel_type
use reference_element_mod,   only: W, S, E, N, B

implicit none

private
public :: restrict_masked_w2_kernel_code

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: restrict_masked_w2_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                                       &
    arg_type(GH_FIELD, GH_REAL, GH_INC,   W2,          mesh_arg=GH_COARSE), & ! coarse field
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  ), & ! fine field
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  ), & ! rmultiplicity
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2, mesh_arg=GH_FINE  )  & ! fine mesh
    /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: restrict_masked_w2_kernel_code
end type restrict_masked_w2_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

contains

  !> @brief Subroutine to perform the W2 restriction operation
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Coarse grid W2 field to compute
  !> @param[in]     fine_field               Fine grid  W2 field to restrict
  !> @param[in]     rmultiplicity            A fine grid W2 field containing the
  !!                                         reciprocal multiplicity of nodes
  !> @param[in]     mask_fine                W2 mask on the fine grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf                      Num of DoFs per cell on both grids
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  subroutine restrict_masked_w2_kernel_code(                  &
                                     nlayers,                 &
                                     cell_map,                &
                                     ncell_fine_per_coarse_x, &
                                     ncell_fine_per_coarse_y, &
                                     ncell_fine,              &
                                     coarse_field,            &
                                     fine_field,              &
                                     rmultiplicity,           &
                                     mask_fine,               &
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
    real(kind=r_def), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_def), intent(in)    :: fine_field(undf_fine)
    real(kind=r_def), intent(in)    :: rmultiplicity(undf_fine)
    real(kind=r_def), intent(in)    :: mask_fine(undf_fine)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ndf, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, ncell_fine_per_coarse_y)

    ! Internal variables
    integer(kind=i_def) :: df, k, x_idx, y_idx, face
    real(kind=r_def)    :: new_coarse
    real(kind=r_def)    :: non_zero_cells

    integer(kind=i_def), parameter :: n_faces = 5
    integer(kind=i_def), parameter :: face_order(n_faces) = [W,S,E,N,B]
    integer(kind=i_def)            :: x_idx_start(n_faces)
    integer(kind=i_def)            :: x_idx_end(n_faces)
    integer(kind=i_def)            :: y_idx_start(n_faces)
    integer(kind=i_def)            :: y_idx_end(n_faces)
    real(kind=r_def)               :: scaling(n_faces)

    !---------------------------------------------------------------------------
    ! Define cells to average over for each df
    !---------------------------------------------------------------------------

    ! The rows and columns forming the cell map match the arrangment
    ! of fine cells within the coarse cell
    !
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

    do face = 1, size(face_order)
      df = face_order(face)

      select case(df)
      case(N)
        ! N edge is first row of cell map
        x_idx_start(df) = 1
        x_idx_end(df)   = ncell_fine_per_coarse_x
        y_idx_start(df) = 1
        y_idx_end(df)   = 1

      case(S)
        ! S edge is last row of cell map
        x_idx_start(df) = 1
        x_idx_end(df)   = ncell_fine_per_coarse_x
        y_idx_start(df) = ncell_fine_per_coarse_y
        y_idx_end(df)   = ncell_fine_per_coarse_y

      case(W)
        ! W edge is first column of cell map
        x_idx_start(df) = 1
        x_idx_end(df)   = 1
        y_idx_start(df) = 1
        y_idx_end(df)   = ncell_fine_per_coarse_y

      case(E)
        ! E edge is last column of cell map
        x_idx_start(df) = ncell_fine_per_coarse_x
        x_idx_end(df)   = ncell_fine_per_coarse_x
        y_idx_start(df) = 1
        y_idx_end(df)   = ncell_fine_per_coarse_y

      case default
        x_idx_start(df) = 1
        x_idx_end(df)   = ncell_fine_per_coarse_x
        y_idx_start(df) = 1
        y_idx_end(df)   = ncell_fine_per_coarse_y

      end select
    end do

    !---------------------------------------------------------------------------
    ! Calculate scaling based on bottom layer of mask
    !---------------------------------------------------------------------------

    scaling(:) = 0.0_r_def
    do face = 1, size(face_order)
      df = face_order(face)
      non_zero_cells = 0.0_r_def

      do y_idx = y_idx_start(df), y_idx_end(df)
        do x_idx = x_idx_start(df), x_idx_end(df)
          non_zero_cells = non_zero_cells + &
                           mask_fine( map_fine(df,cell_map(x_idx,y_idx)) )
        end do
      end do

      if (non_zero_cells > 0.1_r_def) then
        scaling(df) = 1.0_r_def / non_zero_cells
      else
        scaling(df) = 1.0_r_def
      end if
    end do

    !---------------------------------------------------------------------------
    ! Horizontal components
    !---------------------------------------------------------------------------
    ! Shared dofs, so use rmultiplicity and increment the coarse field.

    do face = 1, 4
      df = face_order(face)

      do k = 0, nlayers-1
        new_coarse = 0.0_r_def

        do y_idx = y_idx_start(df), y_idx_end(df)
          do x_idx = x_idx_start(df), x_idx_end(df)
            new_coarse = new_coarse + ( rmultiplicity(map_fine(df,cell_map(x_idx,y_idx))+k) &
                                        * mask_fine(map_fine(df,cell_map(x_idx,y_idx)))     &
                                        * fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) )
          end do
        end do
        coarse_field(map_coarse(df)+k) = coarse_field(map_coarse(df)+k) + new_coarse * scaling(df)
      end do
    end do

    !---------------------------------------------------------------------------
    ! Vertical components
    !---------------------------------------------------------------------------
    ! Only visit dofs once, so don't use rmultiplicity and set coarse field equal
    ! rather than incrementing.
    !
    ! Only do bottom value of cell
    ! Loop over an extra layer to get the very top
    df = B
    do k = 0, nlayers
      new_coarse = 0.0_r_def

      do y_idx = y_idx_start(df), y_idx_end(df)
        do x_idx = x_idx_start(df), x_idx_end(df)
          new_coarse = new_coarse + ( mask_fine(map_fine(df,cell_map(x_idx,y_idx)))    &
                                  * fine_field(map_fine(df,cell_map(x_idx,y_idx))+k) )
        end do
      end do
      coarse_field(map_coarse(df)+k) = new_coarse * scaling(df)
    end do

  end subroutine restrict_masked_w2_kernel_code

end module restrict_masked_w2_kernel_mod
