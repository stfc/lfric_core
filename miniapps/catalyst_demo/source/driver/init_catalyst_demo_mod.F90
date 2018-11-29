!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief init functionality for demonstrating catalyst

!> @details Handles init of prognostic and coordinate fields

module init_catalyst_demo_mod

  use assign_coordinate_field_mod,    only : assign_coordinate_field
  use constants_mod,                  only : i_def
  use field_mod,                      only : field_type, &
                                             write_diag_interface, &
                                             checkpoint_interface, &
                                             restart_interface
  use finite_element_config_mod,      only : element_order

  use fs_continuity_mod,              only : W0, W2, W3, Wtheta
  use gw_init_fields_alg_mod,         only : gw_init_fields_alg
  use log_mod,                        only : log_event, log_scratch_space, &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR
  use restart_control_mod,            only : restart_type
  use gravity_wave_constants_config_mod,only : b_space,                           &
                                               gravity_wave_constants_b_space_w0, &
                                               gravity_wave_constants_b_space_w3, &
                                               gravity_wave_constants_b_space_wtheta
  use runtime_constants_mod,          only : create_runtime_constants
  use output_config_mod,              only : write_xios_output
  use io_mod,                         only : xios_write_field_face, &
                                             xios_write_field_node, &
                                             checkpoint_xios, &
                                             checkpoint_netcdf, &
                                             restart_netcdf, &
                                             restart_xios

  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use function_space_chain_mod,       only : function_space_chain_type

  use init_multigrid_mesh_mod,        only : mesh_ids
  use multigrid_config_mod,           only : l_multigrid, &
                                             ugrid,       &
                                             multigrid_chain_nitems
  implicit none


  contains

  subroutine init_catalyst_demo( mesh_id, chi, multigrid_function_space_chain, &
                                 wind, pressure, buoyancy, restart )

    integer(i_def),                  intent(in)  :: mesh_id
    type(function_space_chain_type), intent(out) :: multigrid_function_space_chain

    ! Prognostic fields
    type( field_type ), intent(inout)        :: wind, pressure, buoyancy
    ! Coordinate field
    type( field_type ), intent(inout)        :: chi(:)
    type(restart_type), intent(in)           :: restart
    integer(i_def)                           :: buoyancy_space

    integer(i_def) :: i

    procedure(write_diag_interface), pointer :: tmp_write_diag_ptr
    procedure(checkpoint_interface), pointer :: tmp_checkpoint_ptr
    procedure(restart_interface), pointer    :: tmp_restart_ptr

    type(function_space_type),  pointer :: function_space => null()

    call log_event( 'catalyst_demo: Initialising miniapp ...', LOG_LEVEL_INFO )

    !===============================================================================
    ! Now create the function space chain for multigrid
    !===============================================================================
    if (l_multigrid) then

      do i=1, multigrid_chain_nitems

        ! Make sure this function_space is in the collection
        function_space => function_space_collection%get_fs( mesh_ids(i), &
                                                            0,           &
                                                            W3 )

        write( log_scratch_space,"(A,I0,A)")                       &
             'Adding function_space id ', function_space%get_id(), &
             ' to multigrid function_space chain'
        call log_event( log_scratch_space, LOG_LEVEL_INFO )

        call multigrid_function_space_chain%add( function_space )
 
      end do
    end if


    ! Create prognostic fields
    select case(b_space)
      case(gravity_wave_constants_b_space_w0)
        buoyancy_space = W0
        call log_event( 'catalyst_demo: Using W0 for buoyancy', LOG_LEVEL_INFO )
      case(gravity_wave_constants_b_space_w3)
        buoyancy_space = W3
        call log_event( 'catalyst_demo: Using W3 for buoyancy', LOG_LEVEL_INFO )
      case(gravity_wave_constants_b_space_wtheta)
        buoyancy_space = Wtheta
        call log_event( 'catalyst_demo: Using Wtheta for buoyancy', LOG_LEVEL_INFO )
      case default
        call log_event( 'catalyst_demo: Invalid buoyancy space', LOG_LEVEL_ERROR )
    end select

    wind = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W2))
    buoyancy = field_type( vector_space = &
               function_space_collection%get_fs(mesh_id, element_order, buoyancy_space) )
    pressure = field_type( vector_space = &
                       function_space_collection%get_fs(mesh_id, element_order, W3) )

    ! Set I/O behaviours for diagnostic output

    if (write_xios_output) then

       ! Fields that are output on the XIOS face domain

       tmp_write_diag_ptr => xios_write_field_face

       call wind%set_write_diag_behaviour(tmp_write_diag_ptr)
       call pressure%set_write_diag_behaviour(tmp_write_diag_ptr)

       if (buoyancy_space == W0) then
         tmp_write_diag_ptr => xios_write_field_node
         call buoyancy%set_write_diag_behaviour(tmp_write_diag_ptr)
       else
         tmp_write_diag_ptr => xios_write_field_face
         call buoyancy%set_write_diag_behaviour(tmp_write_diag_ptr)
       end if
       
    end if

    ! Set I/O behaviours for checkpoint / restart

    if (restart%use_xios()) then

      ! Use XIOS for checkpoint / restart

      tmp_checkpoint_ptr => checkpoint_xios
      tmp_restart_ptr => restart_xios

      call log_event( 'catalyst_demo: Using XIOS for checkpointing...', LOG_LEVEL_INFO )

    else

      ! Use old checkpoint and restart methods

      tmp_checkpoint_ptr => checkpoint_netcdf
      tmp_restart_ptr => restart_netcdf

      call log_event( 'catalyst_demo: Using NetCDF for checkpointing...', LOG_LEVEL_INFO )

    end if

    call wind%set_checkpoint_behaviour(tmp_checkpoint_ptr)
    call pressure%set_checkpoint_behaviour(tmp_checkpoint_ptr)
    call buoyancy%set_checkpoint_behaviour(tmp_checkpoint_ptr)

    call wind%set_restart_behaviour(tmp_restart_ptr)
    call pressure%set_restart_behaviour(tmp_restart_ptr)
    call buoyancy%set_restart_behaviour(tmp_restart_ptr)

    ! Create runtime_constants object. This in turn creates various things
    ! needed by the timestepping algorithms such as mass matrix operators, mass
    ! matrix diagonal fields and the geopotential field
    call create_runtime_constants(mesh_id, chi)

    ! Initialise prognostic fields
    call gw_init_fields_alg(wind, pressure, buoyancy, restart)

    nullify( tmp_write_diag_ptr, tmp_checkpoint_ptr, &
             tmp_restart_ptr, function_space )

    call log_event( 'catalyst_demo: Miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_catalyst_demo

end module init_catalyst_demo_mod
