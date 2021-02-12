!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which computes the departure distances for cell faces in the
!>        vertical direction.
!> @details The Cosmic scheme updates density in the x, y and z directions
!>          separately.
!>          This code calculates the distance which is swept through a cell
!>          in the z direction during one timestep. The arrival point is the
!>          cell face and the departure point is calculated. Options for the A single Euler
!>          calculation of the departure point are a single Euler timestep, the
!>          midpoint rule or trapezoidal rule. This kernel returns the distance
!>          between the arrival and departure point for each cell face in the
!>          vertical and this value is positive if the w wind (radial wind) is
!>          positive, i.e. increasing in height.

module vertical_trapezoidal_kernel_mod

use argument_mod,  only : arg_type, func_type,                  &
                          GH_FIELD, GH_INC, GH_READ,            &
                          GH_BASIS, CELLS
use fs_continuity_mod, only : W0, W2, W3
use constants_mod, only : r_def, i_def
use kernel_mod,    only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vertical_trapezoidal_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                  &
       arg_type(GH_FIELD,   GH_INC, W2),                               &
       arg_type(GH_FIELD,   GH_READ,  W2),                             &
       arg_type(GH_FIELD,   GH_READ,  W2)                              &
       /)
  integer(kind=i_def) :: iterates_over = CELLS
contains
  procedure, nopass :: vertical_trapezoidal_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public vertical_trapezoidal_code
contains

!> @brief Kernel which computes the departure distances for cell faces in the
!>        vertical direction.
!! @param[in]  nlayers             The number of layers
!! @param[in,out]  dep_pts_z       The departure distances in the vertical
!! @param[in]  u_n                 The wind field at time level n
!! @param[in]  u_np1               The wind field at time level n+1
!! @param[in]  ndf_w2              The number of degrees of freedom per cell
!! @param[in]  undf_w2             The number of unique degrees of freedom
!! @param[in]  map_w2              The dofmap for the cell at the base of the column
subroutine vertical_trapezoidal_code(  nlayers,              &
                                       dep_pts_z,            &
                                       u_n,                  &
                                       u_np1,                &
                                       ndf_w2,               &
                                       undf_w2,              &
                                       map_w2 )

  use biperiodic_deppts_mod,       only : calc_dep_point,           &
                                          calc_vertical_trapezoidal
  use biperiodic_deppt_config_mod, only : n_dep_pt_iterations,      &
                                          vertical_method
  use timestepping_config_mod,     only : dt

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                       :: nlayers
  integer(kind=i_def), intent(in)                       :: ndf_w2
  integer(kind=i_def), intent(in)                       :: undf_w2
  integer(kind=i_def), dimension(ndf_w2), intent(in)    :: map_w2
  real(kind=r_def), dimension(undf_w2), intent(in)      :: u_n
  real(kind=r_def), dimension(undf_w2), intent(in)      :: u_np1
  real(kind=r_def), dimension(undf_w2), intent(inout)   :: dep_pts_z

  integer(kind=i_def) :: k

  integer(kind=i_def) :: nCellEdges
  real(kind=r_def)    :: xArrival
  real(kind=r_def),allocatable :: u_n_local(:)
  real(kind=r_def),allocatable :: u_np1_local(:)

  xArrival = 1.0_r_def

  nCellEdges = nlayers+1
  allocate(u_n_local(1:nCellEdges))
  allocate(u_np1_local(1:nCellEdges))

  u_n_local   = -99.0_r_def
  u_np1_local = -99.0_r_def


  ! Extract and fill arrays u_n_local and u_np1_local from global variables
  ! u_n and u_np1.
  do k=1,nlayers-1
    u_n_local(k+1) = u_n(map_w2(5)+k)
    u_np1_local(k+1) = u_np1(map_w2(5)+k)
  end do
  ! Apply vertical boundary conditions.
  u_n_local(1) = 0.0_r_def
  u_np1_local(1) = 0.0_r_def
  u_n_local(nCellEdges) = 0.0_r_def
  u_np1_local(nCellEdges) = 0.0_r_def

  ! Apply vertical boundary conditions to the departure points.
  dep_pts_z( map_w2(5) ) =  0.0_r_def
  dep_pts_z( map_w2(6)+nlayers-1 ) =  0.0_r_def

  ! Loop over all layers except the bottom layer.
  ! This code is hard-wired to work with 6 W2 dofs per cell where dof=5 is the
  ! vertical dof at the bottom of the cell. This code forms part of Cosmic which
  ! is designed to work only for 6 W2 dofs per cell.
  do k=1,nlayers-1
    xArrival = real(k,r_def)
    dep_pts_z( map_w2(5) + k ) =  calc_vertical_trapezoidal( xArrival,             &
                                                             nCellEdges,           &
                                                             u_n_local,            &
                                                             u_np1_local,          &
                                                             dt,                   &
                                                             vertical_method,      &
                                                             n_dep_pt_iterations )
  end do

  deallocate(u_n_local)
  deallocate(u_np1_local)

end subroutine vertical_trapezoidal_code

end module vertical_trapezoidal_kernel_mod
