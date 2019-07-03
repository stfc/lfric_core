!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Controls grid related information used by the model

module gungho_grid_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use global_mesh_collection_mod,     only : global_mesh_collection, &
                                             global_mesh_collection_type
  use init_mesh_mod,                  only : init_mesh
  use init_fem_mod,                   only : init_fem
  use mpi_mod,                        only : get_comm_size, get_comm_rank

  implicit none

  private
  public initialise_grid

contains

  !> @brief Initialises grid related information used by the model
  !> @param [inout] mesh_id The identifier of the primary mesh
  !> @param [inout] twod_mesh_id The identifier of the primary 2d mesh
  !> @param [inout] chi A size 3 array of fields holding the coordinates of the mesh
  subroutine initialise_grid(mesh_id, twod_mesh_id, chi)

    implicit none

    integer(i_def), intent(inout) :: mesh_id, twod_mesh_id
    type(field_type), intent(inout) :: chi(3)

    integer(i_def) :: total_ranks, local_rank
    integer(i_def) :: i

    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    ! Get the rank information from the virtual machine
    total_ranks = get_comm_size()
    local_rank  = get_comm_rank()

    ! Create the mesh
    call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

    ! Create FEM specifics (function spaces and chi field)
    call init_fem(mesh_id, chi)

    ! Full global meshes no longer required, so reclaim
    ! the memory from global_mesh_collection
    deallocate(global_mesh_collection)

  end subroutine initialise_grid

end module gungho_grid_mod
