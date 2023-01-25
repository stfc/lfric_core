!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Remap a scalar field (W3 or Wtheta) from the standard cubed sphere mesh
!!        to the extended cubed sphere
module remap_on_extended_mesh_kernel_mod

use kernel_mod,        only: kernel_type
use argument_mod,      only: arg_type, func_type,       &
                             GH_FIELD, GH_SCALAR,       &
                             GH_REAL, GH_LOGICAL,       &
                             GH_READ, GH_WRITE,         &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             ANY_DISCONTINUOUS_SPACE_3, &
                             GH_BASIS, CELL_COLUMN,     &
                             STENCIL, CROSS2D,          &
                             GH_EVALUATOR
use constants_mod,     only: r_def, i_def, l_def, LARGE_REAL_POSITIVE
use fs_continuity_mod, only: Wchi
use log_mod,           only: log_event, &
                             LOG_LEVEL_ERROR

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: remap_on_extended_mesh_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                                           &
       arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                   &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)), &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi),                                        &
       arg_type(GH_FIELD*3, GH_REAL,    GH_READ,  Wchi, STENCIL(CROSS2D)),                      &
       arg_type(GH_FIELD,   GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3, STENCIL(CROSS2D)), &
       arg_type(GH_SCALAR,  GH_LOGICAL, GH_READ)                                                &
       /)
  type(func_type) :: meta_funcs(1) = (/ &
       func_type(Wchi, GH_BASIS)        &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: remap_on_extended_mesh_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: remap_on_extended_mesh_code
contains

!> @brief Remap a W3 or Wtheta scalar field from the standard cubed sphere mesh
!!        to the extended cubed sphere mesh.
!> @param[in]     nlayers          Number of layers
!> @param[in,out] remap_field      Field to compute remap values in
!> @param[in]     field            Field to map from (on original mesh)
!> @param[in]     stencil_size     Size of the x-stencil (number of cells)
!> @param[in]     stencil          Dofmaps for the x-stencil
!> @param[in]     max_length       Maximum stencil branch length
!> @param[in]     alpha_ext        alpha coordinate on extended mesh
!> @param[in]     beta_ext         beta coordinate on extended mesh
!> @param[in]     height_ext       height coordinate on extended mesh
!> @param[in]     alpha            alpha coordinate on original mesh
!> @param[in]     beta             beta coordinate on original mesh
!> @param[in]     height           height coordinate on original mesh
!> @param[in]     wx_stencil_size  Size of the Wx-stencil (number of cells)
!> @param[in]     wx_stencil       Dofmaps for the Wx-stencil
!> @param[in]     wx_max_length    Maximum stencil branch length
!> @param[in]     panel_id         Id of cubed sphere panel
!> @param[in]     pid_stencil_size Size of the panelid-stencil (number of cells)
!> @param[in]     pid_stencil      Dofmaps for the panelid-stencil
!> @param[in]     pid_max_length   Maximum stencil branch length
!> @param[in]     linear_remap     Use a linear (=true) or cubic (=false) remapping
!> @param[in]     ndf_ws           Number of degrees of freedom per cell for scalar fields
!> @param[in]     undf_ws          Number of unique degrees of freedom for scalar fields
!> @param[in]     map_ws           Dofmap for the cell at the base of the column for scalar fields
!> @param[in]     ndf_wx           Number of degrees of freedom per cell for coord fields
!> @param[in]     undf_wx          Number of unique degrees of freedom for coord fields
!> @param[in]     map_wx           Dofmap for the cell at the base of the column for coord fields
!! @param[in]     basis_wx         Basis functions of coordinate field evaluated at nodes of scalar fields
!> @param[in]     ndf_pid          Number of degrees of freedom per cell for the panel id field
!> @param[in]     undf_pid         Number of unique degrees of freedom for the panel id field
!> @param[in]     map_pid          Dofmap for the cell at the base of the column for the panel id field
subroutine remap_on_extended_mesh_code(nlayers,                                       &
                                       remap_field,                                   &
                                       field,                                         &
                                       stencil_size, stencil, max_length,             &
                                       alpha_ext, beta_ext, height_ext,               &
                                       alpha, beta, height,                           &
                                       wx_stencil_size, wx_stencil, wx_max_length,    &
                                       panel_id,                                      &
                                       pid_stencil_size, pid_stencil, pid_max_length, &
                                       linear_remap,                                  &
                                       ndf_ws, undf_ws, map_ws,                       &
                                       ndf_wx, undf_wx, map_wx, basis_wx,             &
                                       ndf_pid, undf_pid, map_pid                     &
                                      )

  use coord_transform_mod, only: alphabetar2xyz, &
                                 xyz2alphabetar

  implicit none

  ! Arguments
  integer(kind=i_def),                     intent(in) :: nlayers
  integer(kind=i_def), dimension(4),       intent(in) :: stencil_size
  integer(kind=i_def), dimension(4),       intent(in) :: wx_stencil_size
  integer(kind=i_def), dimension(4),       intent(in) :: pid_stencil_size
  integer(kind=i_def),                     intent(in) :: max_length
  integer(kind=i_def),                     intent(in) :: wx_max_length
  integer(kind=i_def),                     intent(in) :: pid_max_length
  integer(kind=i_def),                     intent(in) :: ndf_ws, undf_ws
  integer(kind=i_def),                     intent(in) :: ndf_wx, undf_wx
  integer(kind=i_def),                     intent(in) :: ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_ws),  intent(in) :: map_ws
  integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  logical(kind=l_def), intent(in) :: linear_remap

  integer(kind=i_def), dimension(ndf_ws,  max_length,    4), intent(in) :: stencil
  integer(kind=i_def), dimension(ndf_wx,  wx_max_length, 4), intent(in) :: wx_stencil
  integer(kind=i_def), dimension(ndf_pid, pid_max_length,4), intent(in) :: pid_stencil

  real(kind=r_def), dimension(1,ndf_wx,ndf_ws), intent(in) :: basis_wx

  real(kind=r_def), dimension(undf_ws),  intent(inout) :: remap_field
  real(kind=r_def), dimension(undf_ws),  intent(in)    :: field
  real(kind=r_def), dimension(undf_wx),  intent(in)    :: alpha_ext, beta_ext, height_ext
  real(kind=r_def), dimension(undf_wx),  intent(in)    :: alpha, beta, height
  real(kind=r_def), dimension(undf_pid), intent(in)    :: panel_id

  integer(kind=i_def)               :: owned_panel, halo_panel, panel, &
                                       panel_edge, df, k,              &
                                       id1, id2, id3, id4, n, nl
  integer(kind=i_def)               :: interp_dir
  integer(kind=i_def), parameter    :: interp_dir_alpha = 1
  integer(kind=i_def), parameter    :: interp_dir_beta  = 2
  integer(kind=i_def), dimension(2) :: ncells_in_stencil

  integer(kind=i_def), dimension(2*max_length-1, 2)         :: stencil_1d
  integer(kind=i_def), dimension(ndf_wx, 2*max_length-1, 2) :: wx_stencil_1d
  integer(kind=i_def), dimension(2*max_length-1, 2)         :: pid_stencil_1d

  real(kind=r_def), dimension(3) :: abh
  real(kind=r_def)               :: x0, x, y, z, w1, w2, w3, w4
  real(kind=r_def), allocatable  :: x1(:), dx(:)
  real(kind=r_def), parameter    :: unit_radius = 1.0_r_def
  integer(kind=i_def)            :: alpha_dir_id, beta_dir_id

  ! Assume the first entry in panel id corresponds to an owned (not halo) cell
  owned_panel = int(panel_id(1),i_def)
  ! Panel id for this column
  halo_panel = int(panel_id(map_pid(1)),i_def)

  ! nl = nlayers for Wtheta fields (ndf_ws=2)
  ! nl = nlayers-1 for W3 fields (ndf_ws=1)
  nl = (nlayers - 1) + (ndf_ws - 1)

  ! Only need to remap if the halo is on a different panel to the
  ! owned cell (otherwise the remaped field is identical to the original
  ! field and should have been picked up by a copy_field call)
  if ( halo_panel /= owned_panel ) then

    ! Compute interpolation direction (alpha or beta)
    ! the first digit of panel_edge is the panel we are interpolating to and the
    ! second is the one we are interpolating from, i.e. panel_edge = 36 means we
    ! are interpolating from panel 6 to panel 3
    panel_edge = 10*owned_panel + halo_panel
    select case (panel_edge)
      case (15, 16, 25, 26, 32, 34, 41, 43, 51, 53, 62, 64)
        interp_dir = interp_dir_alpha
      case (12, 14, 21, 23, 35, 36, 45, 46, 52, 54, 61, 63)
        interp_dir = interp_dir_beta
      case default
        call log_event('Invalid panel edge',LOG_LEVEL_ERROR)
    end select

    ! The stencils are ordered (W,S,E,N) on their local panel with (W,E) in the
    ! alpha direction and (S,N) in the beta direction. However, this doesn't
    ! necessarily correspond to the same direction on the owned panel so we
    ! need to potentially rotate the directions.
    ! i.e on the halo of panel 1 corresponding to panel 4 we want to
    ! interpolate in the beta direction of panel 1 but this corresponds to
    ! alpha direction on panel 4 so we set alpha_dir_id = 2 here
    select case( panel_edge )
      case (14, 16, 23, 25, 32, 35, 41, 46, 52, 53, 61, 64)
        ! Change in orientation of alpha and beta directions
        alpha_dir_id = 2
        beta_dir_id  = 1
      case default
        ! Default values for when there is no change in orientation
        alpha_dir_id = 1
        beta_dir_id  = 2
    end select

    ! Compute combined 1D stencils
    ! for a Cross stencil of the form:
    !      | 7 |
    !      | 5 |
    !  | 2 | 1 | 4 | 6 |
    !      | 3 |
    ! Stored as stencil =[
    ! 1, 2;
    ! 1, 3;
    ! 1, 4, 6;
    ! 1, 5, 7]
    ! Then the new stencil_1d = [
    ! 1, 2, 4, 6;
    ! 1, 3, 5, 7]

    ncells_in_stencil(1) = stencil_size(alpha_dir_id) + stencil_size(alpha_dir_id+2) - 1
    ncells_in_stencil(2) = stencil_size(beta_dir_id)  + stencil_size(beta_dir_id+2)  - 1

    stencil_1d = 0
    do n = 1, stencil_size(alpha_dir_id)
      stencil_1d(n,1) = stencil(1,n,alpha_dir_id)
    end do
    do n = 1, stencil_size(alpha_dir_id+2)-1
      stencil_1d(n+stencil_size(alpha_dir_id),1) = stencil(1,n+1,alpha_dir_id+2)
    end do
    do n = 1, stencil_size(beta_dir_id)
      stencil_1d(n,2) = stencil(1,n,beta_dir_id)
    end do
    do n = 1, stencil_size(beta_dir_id+2)-1
      stencil_1d(n+stencil_size(beta_dir_id),2) = stencil(1,n+1,beta_dir_id+2)
    end do

    wx_stencil_1d = 0
    do n = 1, wx_stencil_size(alpha_dir_id)
      wx_stencil_1d(:,n,1) = wx_stencil(:,n,alpha_dir_id)
    end do
    do n = 1, wx_stencil_size(alpha_dir_id+2)-1
      wx_stencil_1d(:,n+wx_stencil_size(alpha_dir_id),1) = wx_stencil(:,n+1,alpha_dir_id+2)
    end do
    do n = 1, wx_stencil_size(beta_dir_id)
      wx_stencil_1d(:,n,2) = wx_stencil(:,n,beta_dir_id)
    end do
    do n = 1, wx_stencil_size(beta_dir_id+2)-1
      wx_stencil_1d(:,n+wx_stencil_size(beta_dir_id),2) = wx_stencil(:,n+1,beta_dir_id+2)
    end do

    pid_stencil_1d = 0
    do n = 1, pid_stencil_size(alpha_dir_id)
      pid_stencil_1d(n,1) = pid_stencil(1,n,alpha_dir_id)
    end do
    do n = 1, pid_stencil_size(alpha_dir_id+2)-1
      pid_stencil_1d(n+pid_stencil_size(alpha_dir_id),1) = pid_stencil(1,n+1,alpha_dir_id+2)
    end do
    do n = 1, pid_stencil_size(beta_dir_id)
      pid_stencil_1d(n,2) = pid_stencil(1,n,beta_dir_id)
    end do
    do n = 1, pid_stencil_size(beta_dir_id+2)-1
      pid_stencil_1d(n+pid_stencil_size(beta_dir_id),2) = pid_stencil(1,n+1,beta_dir_id+2)
    end do

    ! Compute interpolation point (x0) in centre of remapped cell
    abh = 0.0_r_def
    do df = 1, ndf_wx
      abh(1) = abh(1) + alpha_ext(map_wx(df))*basis_wx(1,df,1)
      abh(2) = abh(2) + beta_ext(map_wx(df))*basis_wx(1,df,1)
    end do
    x0 = abh(interp_dir)

    ! Now compute all points on original mesh in coordinates of extended mesh
    allocate( x1(ncells_in_stencil(interp_dir)), &
              dx(ncells_in_stencil(interp_dir)) )

    do n = 1,ncells_in_stencil(interp_dir)

      panel = int(panel_id(pid_stencil_1d(n, interp_dir)), i_def)

      if ( halo_panel /= panel ) then
        ! If point is on a different panel to the halo panel (i.e we have gone
        ! 'around' the cubed sphere corner) then discard that point (by setting the distance to
        ! some large number
        x1(n) = LARGE_REAL_POSITIVE
      else
        ! Otherwise compute the coordinates of the point in terms of the
        ! owned_panel
        abh = 0.0_r_def
        do df = 1, ndf_wx
          abh(1) = abh(1) + alpha(wx_stencil_1d(df, n, interp_dir))*basis_wx(1,df,1)
          abh(2) = abh(2) + beta(wx_stencil_1d(df, n, interp_dir))*basis_wx(1,df,1)
          abh(3) = abh(3) + height(wx_stencil_1d(df, n, interp_dir))*basis_wx(1,df,1)
        end do
        abh(3) = abh(3) + unit_radius
        ! Convert (alpha,beta,r) on panel halo_panel to xyz
        call alphabetar2xyz(abh(1), abh(2), abh(3), &
                            panel, x, y, z)
        ! Now convert back to (alpha, beta, r) on owned panel
        call xyz2alphabetar(x, y, z, owned_panel, &
                            abh(1), abh(2), abh(3) )
        x1(n) = abh(interp_dir)
      end if
    end do

    ! Distance between interpolation point x0 and data points x1
    dx = abs(x1 - x0)
    ! Now find two points (for linear interpolation) in x1 that are closest to x0
    id1 = minloc(dx,1)
    dx(id1) = LARGE_REAL_POSITIVE
    id2 = minloc(dx,1)

    if ( linear_remap ) then
      ! Compute weights for linear interpolation
      w1 = (x0 - x1(id2))/(x1(id1) - x1(id2))
      w2 = (x0 - x1(id1))/(x1(id2) - x1(id1))
      ! Remap all dofs in the column using the same weights
      do k = 0, nl
        remap_field(map_ws(1)+k) = w1 * field( stencil_1d(id1,interp_dir)+k) &
                                 + w2 * field( stencil_1d(id2,interp_dir)+k)
      end do
    else
      ! Compute weights for cubic interpolation
      ! First find the next two closest points in stencil
      dx(id2) = LARGE_REAL_POSITIVE
      id3 = minloc(dx,1)
      dx(id3) = LARGE_REAL_POSITIVE
      id4 = minloc(dx,1)

      ! Now compute the weights
      w1 = (x0 - x1(id2))/(x1(id1) - x1(id2)) * (x0 - x1(id3))/(x1(id1) - x1(id3)) * (x0 - x1(id4))/(x1(id1) - x1(id4))
      w2 = (x0 - x1(id1))/(x1(id2) - x1(id1)) * (x0 - x1(id3))/(x1(id2) - x1(id3)) * (x0 - x1(id4))/(x1(id2) - x1(id4))
      w3 = (x0 - x1(id1))/(x1(id3) - x1(id1)) * (x0 - x1(id2))/(x1(id3) - x1(id2)) * (x0 - x1(id4))/(x1(id3) - x1(id4))
      w4 = (x0 - x1(id1))/(x1(id4) - x1(id1)) * (x0 - x1(id2))/(x1(id4) - x1(id2)) * (x0 - x1(id3))/(x1(id4) - x1(id3))
      ! Remap all dofs in the column using the same weights
      do k = 0, nl
        remap_field(map_ws(1)+k) = w1 * field( stencil_1d(id1,interp_dir)+k) &
                                 + w2 * field( stencil_1d(id2,interp_dir)+k) &
                                 + w3 * field( stencil_1d(id3,interp_dir)+k) &
                                 + w4 * field( stencil_1d(id4,interp_dir)+k)

      end do
    end if

    deallocate( x1, dx )

  else
    ! Halo value is on same panel as owned value so just copy the field
    do k = 0, nl
      remap_field(map_ws(1)+k) = field(map_ws(1)+k)
    end do

  end if

end subroutine remap_on_extended_mesh_code

end module remap_on_extended_mesh_kernel_mod
