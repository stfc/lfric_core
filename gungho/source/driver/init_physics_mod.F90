!-------------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
!-------------------------------------------------------------------------------
!> @brief init functionality for physics

!> @details Handles initialization of prognostic fields

module init_physics_mod

  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta
  use init_prognostic_fields_alg_mod, only : init_prognostic_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,         &
                                             LOG_LEVEL_ERROR
  use restart_control_mod,            only : restart_type
  use formulation_config_mod,         only : transport_only
  use transport_config_mod,           only : scheme, &
                                             operators, &
                                             transport_scheme_method_of_lines, &
                                             transport_operators_fv
  use mr_indices_mod,                 only : nummr
  use runtime_constants_mod,          only : create_runtime_constants
  use map_physics_fields_alg_mod,     only : map_physics_fields_alg

  implicit none


  contains

  subroutine init_physics(mesh_id, u, rho, theta, rho_in_wth,         &
                          u1_in_w3, u2_in_w3, u3_in_w3, theta_in_w3,  &
                          p_in_w3, p_in_wth)

    integer(i_def), intent(in)               :: mesh_id
    ! prognostic fields
    type( field_type ), intent(inout)        :: u, rho, theta, rho_in_wth
    type( field_type ), intent(inout)        :: u1_in_w3, u2_in_w3, u3_in_w3, theta_in_w3
    type( field_type ), intent(inout)        :: p_in_w3, p_in_wth

    integer(i_def) :: theta_space
    call log_event( 'Physics: initialisation...', LOG_LEVEL_INFO )
    

    theta_space=theta%which_function_space()
    
    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta', LOG_LEVEL_ERROR )
    end if


    ! Note that the fields generated here should realy be finite volume 
    ! or finite difference fields, but we exploit the lowest order W3/Wtheta finite
    ! element construction.  This may break without rehabilitation, so this 
    ! need doing properly.  Ticket #1044 is opened to address this.
    u1_in_w3 = field_type( vector_space = & 
       function_space_collection%get_fs(mesh_id, element_order, W3))

    u2_in_w3 = field_type( vector_space = & 
       function_space_collection%get_fs(mesh_id, element_order, W3))

    u3_in_w3 = field_type( vector_space = & 
       function_space_collection%get_fs(mesh_id, element_order, W3))

    theta_in_w3 = field_type( vector_space = & 
       function_space_collection%get_fs(mesh_id, element_order, W3))

    p_in_w3 = field_type( vector_space = & 
       function_space_collection%get_fs(mesh_id, element_order, W3))

    p_in_wth = field_type( vector_space = & 
       function_space_collection%get_fs(mesh_id, element_order, theta_space))

    call map_physics_fields_alg( u, rho, theta, rho_in_wth,                    &
                                 u1_in_w3, u2_in_w3, u3_in_w3, theta_in_w3,    &
                                 p_in_w3, p_in_wth)

    call log_event( 'Physics initialised', LOG_LEVEL_INFO )

  end subroutine init_physics

end module init_physics_mod
