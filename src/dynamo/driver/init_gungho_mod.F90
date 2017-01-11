!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown,
! Met Office and NERC 2014.
! However, it has been created with the help of the GungHo Consortium,
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief init functionality for gungho

!> @details Handles mesh creation and function space creation

module init_gungho_mod

  use constants_mod,                  only : i_def
  use ESMF
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use log_mod,                        only : log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO
  use global_mesh_collection_mod,     only : global_mesh_collection_type, &
                                             global_mesh_collection
  use mesh_collection_mod,            only : mesh_collection_type, &
                                             mesh_collection
  use set_up_mod,                     only : set_up


  implicit none


  contains

  subroutine init_gungho(mesh_id, local_rank, total_ranks)

    integer(i_def), intent(out)   :: mesh_id
    integer(i_def), intent(in)    :: total_ranks
    integer(i_def), intent(in)    :: local_rank

    ! Create top level collections
    allocate( global_mesh_collection, &
              source = global_mesh_collection_type() )

    allocate( mesh_collection, &
              source = mesh_collection_type() )

    allocate( function_space_collection, &
              source = function_space_collection_type() )

    ! Set up mesh
    call set_up(local_rank, total_ranks, mesh_id)

    ! Full global meshes no longer required, so reclaim the memory
    ! from global_mesh_collection
    write(log_scratch_space,'(A)') &
        "Gungho: Meshes created, purging global mesh collection."
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

    deallocate(global_mesh_collection)

    write(log_scratch_space,'(A,I0)') "Gungho: Partition mesh, id=", mesh_id
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

    call log_event( 'Gungho initialised', LOG_LEVEL_INFO )

  end subroutine init_gungho

end module init_gungho_mod
