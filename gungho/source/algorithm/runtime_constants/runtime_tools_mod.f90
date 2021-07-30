!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Provides some generic routines used in setting up runtime constants
!>
!> @details This module provides a few tools common to each of the sets of
!>          runtime constants groups:
!>          - Global mesh ID lists are provided here as they are needed
!>            throughout run-time by the runtime constants getter functions.
!>          - Meshes are given labels to describe what type of mesh they are.
!>            This is used to determine which constants to set up on which mesh,
!>            as many constants aren't needed on every mesh. The enumerated
!>            labels are provided here.
!>          Having these in a separate module is necessary for avoiding circular
!>          dependency trees and duplicated code.
module runtime_tools_mod

  use constants_mod,  only: i_def
  use log_mod,        only: log_event, LOG_LEVEL_ERROR

  implicit none

  private

  ! Mesh IDs
  integer(kind=i_def), allocatable :: global_mesh_id_list(:)

  ! Enumerated types to label what type the mesh is (randomly selected integers)
  integer(kind=i_def), parameter, public :: primary_mesh_label = 19
  integer(kind=i_def), parameter, public :: shifted_mesh_label = 964
  integer(kind=i_def), parameter, public :: double_level_mesh_label = 209
  integer(kind=i_def), parameter, public :: multigrid_mesh_label = 624
  integer(kind=i_def), parameter, public :: twod_mesh_label = 71
  integer(kind=i_def), parameter, public :: extra_mesh_label = 117

  ! Public functions to create and access the module contents
  public :: init_mesh_id_list
  public :: final_mesh_id_list
  public :: find_mesh_index

contains
  !> @brief Subroutine to initialise mesh ID list
  !> @param[in] mesh_id_list          List of mesh IDs
  subroutine init_mesh_id_list(mesh_id_list)

    implicit none

    integer(kind=i_def), intent(in) :: mesh_id_list(:)

    ! Internal variables
    integer(kind=i_def) :: i, number_of_meshes

    number_of_meshes = size(mesh_id_list)

    allocate(global_mesh_id_list(number_of_meshes))

    do i = 1, number_of_meshes
      global_mesh_id_list(i) = mesh_id_list(i)
    end do

  end subroutine init_mesh_id_list

  !> @brief Deallocates the mesh ID list
  subroutine final_mesh_id_list()
    implicit none

    if (allocated(global_mesh_id_list)) deallocate(global_mesh_id_list)

  end subroutine final_mesh_id_list

  !> @brief Returns the mesh index in the mesh list from the mesh_id
  !> @param[in] mesh_id
  !> @return The mesh index
  function find_mesh_index(mesh_id) result(mesh_index)
    ! TODO: This routine should be removed in #2790. It should be possible
    ! to return the appropriate constant without relying on a mesh_index

    implicit none
    integer(kind=i_def), intent(in) :: mesh_id
    integer(kind=i_def)             :: i, mesh_index

    mesh_index = 0_i_def
    do i = 1, size(global_mesh_id_list)
      if ( global_mesh_id_list(i) == mesh_id ) then
        mesh_index = i
        exit
      end if
    end do

    if ( mesh_index == 0_i_def ) then
      call log_event( "mesh_id cannot be found in mesh_id_list", LOG_LEVEL_ERROR )
    end if
  end function find_mesh_index

end module runtime_tools_mod
