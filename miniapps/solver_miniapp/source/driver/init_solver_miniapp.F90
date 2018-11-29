!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief init functionality for the miniapp skeleton

!> @details Handles init of prognostic fields and through the call to 
!>          runtime_csontants the coordinate fields and fem operators

module init_solver_miniapp_mod

  use constants_mod,                  only : i_def, r_def
  use field_mod,                      only : field_type
  use field_vector_mod,               only : field_vector_type
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use fs_continuity_mod,              only : W0, W3
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR, &
                                             log_scratch_space
  use psykal_lite_mod,                only : invoke_set_field_scalar
  use runtime_constants_mod,          only : create_runtime_constants
  implicit none


contains
  !> initialise the fields and field vector for the miniapp.
  !> @param[in] mesh_id Integer, the id of the mesh
  !> @param[inout] chi An size 3 array of fields holding the coordinates of the mesh
  !> @param[inout] The field vector which is to be initialised.
  subroutine init_solver_miniapp(mesh_id, chi, fv )
    implicit none
    integer(i_def), intent(in)             :: mesh_id
    type( field_type ), intent(inout)      :: chi(:)
    ! prognostic fields
    type( field_vector_type ), intent(inout) :: fv
    type( field_type )                       :: f1, f2

    
    call log_event( 'solver miniapp: initialisation...', LOG_LEVEL_INFO )
    
    
    ! Create prognostic fields 
    ! Create a field in the W0 function space (fully continuous field)
    f1 = field_type( vector_space = &
         function_space_collection%get_fs(mesh_id, element_order, W0) )
    ! Create a field in the W3 function space (fully discontinuous field)
    f2 = field_type( vector_space = &
         function_space_collection%get_fs(mesh_id, element_order, W3) )
    
    ! set the fields to scalars
    call invoke_set_field_scalar(0.5_r_def,f1)
    call invoke_set_field_scalar(1.0_r_def,f2)
    ! make the vector
    fv = field_vector_type(2_i_def)
    call fv%import_field(f1,1)
    call fv%import_field(f2,2)
    
    write(log_scratch_space,'(A,E16.8)') "W0, W3 Fvector initiialised to 0.5/1.0 norm=",fv%norm()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    
    
    call create_runtime_constants(mesh_id, chi)
    
    call log_event( 'solver miniapp initialised', LOG_LEVEL_INFO )
    
  end subroutine init_solver_miniapp
  
end module init_solver_miniapp_mod
