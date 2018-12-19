!-------------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Init functionality for physics
!> @details Handles initialization of prognostic fields
module init_physics_mod

  use constants_mod,                  only : i_def, l_def
  use field_mod,                      only : field_type, &
                                             write_diag_interface, &
                                             checkpoint_interface, &
                                             restart_interface
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection
  use field_collection_mod,           only : field_collection_type
  use fs_continuity_mod,              only : W0, W1, W2, W3, Wtheta
  use function_space_mod,             only : function_space_type
  use init_cloud_twod_fields_alg_mod, only : init_cloud_twod_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO,         &
                                             LOG_LEVEL_ERROR
  use restart_control_mod,            only : restart_type
  use formulation_config_mod,         only : transport_only, &
                                             use_moisture


  use transport_config_mod,           only : scheme, &
                                             operators, &
                                             transport_scheme_method_of_lines, &
                                             transport_operators_fv
  use mr_indices_mod,                 only : nummr
  use runtime_constants_mod,          only : create_runtime_constants

  implicit none

contains
  !>@brief Routine to initialise the field objects required by the physics
  !> @param[in] mesh_id Identifier of the mesh
  !> @param[in] twod_mesh_id Identifier of the 2D (surface) mesh
  !> @param[in] restart checkpoint/restart type with timestepping information 
  !> @param[in,out] u Wind field
  !> @param[in,out] exner Exner pressure field
  !> @param[in,out] rho Density field
  !> @param[in,out] theta Potential temperature field
  !> @param[out]   derived_fields Collection of FD fields derived from FE fields 
  !> @param[out]   cloud_fields Collection of FD cloud fields
  !> @param[out]   twod_fields Collection of two fields
  !> @param[out]   physics_incs Collection of physics increments
  subroutine init_physics(mesh_id, twod_mesh_id, restart,             &
                          u, exner, rho, theta,                       &
                          derived_fields, cloud_fields, twod_fields,  &
                          physics_incs)

    implicit none

    integer(i_def), intent(in)               :: mesh_id
    type(restart_type), intent(in)           :: restart
    integer(i_def), intent(in)               :: twod_mesh_id
    ! Prognostic fields
    type( field_type ), intent(inout)        :: u, exner, rho, theta

    ! Collections of fields
    type(field_collection_type), intent(out) :: twod_fields
    type(field_collection_type), intent(out) :: cloud_fields
    type(field_collection_type), intent(out) :: derived_fields
    type(field_collection_type), intent(out) :: physics_incs

    ! pointers to vector spaces
    type(function_space_type), pointer  :: vector_space => null()

    integer(i_def) :: theta_space
    logical(l_def) :: checkpoint_restart_flag

    call log_event( 'Physics: initialisation...', LOG_LEVEL_INFO )
    
    theta_space=theta%which_function_space()
    
    if (theta_space /= Wtheta)then
      call log_event( 'Physics: requires theta variable to be in Wtheta', LOG_LEVEL_ERROR )
    end if

    if (element_order > 0)then
      call log_event( 'Physics: requires lowest order elements', LOG_LEVEL_ERROR )
    end if

    !========================================================================
    ! Here we create some field collections
    !========================================================================
    derived_fields  =  field_collection_type(name='derived_fields')
    checkpoint_restart_flag = .false. 
   
    ! Wtheta fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)

    call add_physics_field(derived_fields, 'w_physics',      vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'rho_in_wth',     vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'wetrho_in_wth',  vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'exner_in_wth',   vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'w_physics_star', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'theta_star',     vector_space, checkpoint_restart_flag, restart)

    ! W3 fields
    vector_space=> function_space_collection%get_fs(mesh_id, 0, W3) 

    call add_physics_field(derived_fields, 'u1_in_w3',      vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u2_in_w3',      vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u3_in_w3',      vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'theta_in_w3',   vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'wetrho_in_w3',  vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u1_in_w3_star', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u2_in_w3_star', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u3_in_w3_star', vector_space, checkpoint_restart_flag, restart)

    ! W2 fields
    vector_space=>function_space_collection%get_fs(mesh_id, 0, W2)

    call add_physics_field(derived_fields, 'u_physics',      vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u_physics_star', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(derived_fields, 'u_star',         vector_space, checkpoint_restart_flag, restart)


    !========================================================================
    ! Here we create some 2d fields for the UM physics
    !========================================================================
    twod_fields = field_collection_type(name='twod_fields')
    vector_space=> function_space_collection%get_fs(twod_mesh_id, 0, W3) 
    checkpoint_restart_flag = .true.

    call add_physics_field(twod_fields, 'tstar',   vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(twod_fields, 'zh',      vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(twod_fields, 'z0msea',  vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(twod_fields, 'ntml',    vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(twod_fields, 'cumulus', vector_space, checkpoint_restart_flag, restart)

    checkpoint_restart_flag = .false.
    call add_physics_field(twod_fields, 'cos_zenith_angle',   vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(twod_fields, 'lit_fraction',       vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(twod_fields, 'stellar_irradiance', vector_space, checkpoint_restart_flag, restart)

    !========================================================================
    ! Here we create some cloud fields
    !========================================================================

    cloud_fields = field_collection_type(name='cloud_fields')
    vector_space=>function_space_collection%get_fs(mesh_id, 0, Wtheta)
    checkpoint_restart_flag = .true.

    call add_physics_field(cloud_fields, 'area_fraction',    vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(cloud_fields, 'ice_fraction',     vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(cloud_fields, 'liquid_fraction',  vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(cloud_fields, 'bulk_fraction',    vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(cloud_fields, 'rhcrit',           vector_space, checkpoint_restart_flag, restart)

    ! Initialise cloud fields
    call init_cloud_twod_fields_alg( cloud_fields, twod_fields, restart )

    !========================================================================
    ! Increment values from individual physics parametrizations
    ! either for use by subsequent parametrizations or as diagnostics
    !========================================================================
    physics_incs = field_collection_type(name='physics_incs')
    vector_space => function_space_collection%get_fs(mesh_id, 0, Wtheta)
    checkpoint_restart_flag = .false. ! no need to dump any of these

    call add_physics_field(physics_incs, 'dt_bl', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(physics_incs, 'dmv_bl', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(physics_incs, 'dt_conv', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(physics_incs, 'dmv_conv', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(physics_incs, 'dtl_mphys', vector_space, checkpoint_restart_flag, restart)
    call add_physics_field(physics_incs, 'dmt_mphys', vector_space, checkpoint_restart_flag, restart)

    call log_event( 'Physics initialised', LOG_LEVEL_INFO ) 

  end subroutine init_physics

  !>@brief Add field to field collection and set its write, 
  !>       checkpoint and restart behaviour 
  !> @param[in,out] field_collection Field collection that 'name' will be added to
  !> @param[in]     name             Name of field to be added to collection 
  !> @param[in]     vector_space     Function space of field to set behaviour for
  !> @param[in]     checkpoint_restart_flag  Flag to allow checkpoint and
  !>                                        restart behaviour of field to be set
  !> @param[in]     restart          Checkpoint/restart type with timestepping
  !>                                information 
  subroutine add_physics_field(field_collection, name, vector_space, &
                               checkpoint_restart_flag, restart)

    use output_config_mod,  only : write_xios_output
    use io_mod,             only : xios_write_field_face, &
                                   checkpoint_xios,       &
                                   checkpoint_netcdf,     &
                                   restart_netcdf,        &
                                   restart_xios

    implicit none
    
    character(*), intent(in)                   :: name
    type(field_collection_type), intent(inout) :: field_collection
    type(function_space_type), intent(in)      :: vector_space
    type(restart_type), intent(in)             :: restart
    logical(l_def), intent(in)                 :: checkpoint_restart_flag

    !Local variables
    type(field_type)                           :: new_field

    ! pointers for xios write interface
    procedure(write_diag_interface), pointer   :: write_diag_behaviour => null()
    procedure(checkpoint_interface), pointer   :: checkpoint_behaviour => null() 
    procedure(restart_interface), pointer      :: restart_behaviour => null()

    ! All physics fields currently require output on faces...
    write_diag_behaviour => xios_write_field_face

    if (restart%use_xios()) then
      checkpoint_behaviour => checkpoint_xios
      restart_behaviour => restart_xios
    else
      checkpoint_behaviour => checkpoint_netcdf
      restart_behaviour => restart_netcdf
    endif

    new_field = field_type( vector_space, name=trim(name) )

    if (write_xios_output) then
      call new_field%set_write_diag_behaviour(write_diag_behaviour)
    end if
    if (checkpoint_restart_flag) then
      call new_field%set_checkpoint_behaviour(checkpoint_behaviour)
      call new_field%set_restart_behaviour(restart_behaviour)
    endif

    call field_collection%add_field(new_field)

  end subroutine add_physics_field

end module init_physics_mod
