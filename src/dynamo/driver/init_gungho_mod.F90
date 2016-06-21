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
  use mesh_mod,                       only : mesh_type
  use set_up_mod,                     only : set_up


  implicit none


  contains

  subroutine init_gungho(mesh_id, local_rank, total_ranks,             &
                         function_space_collection)

    integer(i_def), intent(inout)                        :: mesh_id

    type(function_space_collection_type), intent(inout)  :: function_space_collection
    integer(i_def), intent(in)                           :: total_ranks
    integer(i_def), intent(in)                           :: local_rank



    ! Set up mesh
    call set_up(local_rank, total_ranks, mesh_id)

    write(log_scratch_space,'(A,I0)') "Gungho:Partition mesh, id=", mesh_id
    call log_event(log_scratch_space,LOG_LEVEL_INFO)

    ! Create top level function space collection

    function_space_collection = function_space_collection_type()

    call log_event( 'Gungho initialised', LOG_LEVEL_INFO )

 


  end subroutine init_gungho

end module init_gungho_mod
