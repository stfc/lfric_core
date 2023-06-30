!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the multires_coupling miniapp.
!>
!> This a base structure from which further physics-dynamics coupling development can
!> be done. A test algorithm can be called in order to test the scalar mapping and
!> setting up fields on different meshes. If the miniapp is not in test mode, the
!> regular gungho timestep is called on the dynamics mesh.
!>
module multires_coupling_driver_mod

  use multires_coupling_config_mod, &
    only : dynamics_mesh_name,      &
           physics_mesh_name,       &
           multires_coupling_mode,  &
           multires_coupling_mode_test
  use variable_fields_mod,                      only : update_variable_fields
  use gungho_modeldb_mod,                       only : modeldb_type
  use gungho_step_mod,                          only : gungho_step
  use constants_mod,                            only : i_def, i_native, &
                                                       str_def, imdi
  use io_config_mod,                            only : write_diag,           &
                                                       use_xios_io,          &
                                                       diagnostic_frequency, &
                                                       nodal_output_on_w3
  use io_context_mod,                           only : io_context_type
  use mesh_collection_mod,                      only : mesh_collection
  use log_mod,                                  only : log_event,        &
                                                       LOG_LEVEL_ALWAYS, &
                                                       LOG_LEVEL_INFO
  use mesh_mod,                                 only : mesh_type
  use model_clock_mod,                          only : model_clock_type
  use mpi_mod,                                  only : mpi_type
  use coupling_test_alg_mod,                    only : coupling_test_alg
  use xios,                                     only : xios_context_finalize
  use multires_coupling_model_mod,              only : initialise_model,          &
                                                       finalise_model,            &
                                                       initialise_infrastructure, &
                                                       finalise_infrastructure
  use gungho_model_data_mod,                    only : create_model_data,     &
                                                       finalise_model_data,   &
                                                       initialise_model_data
  use multires_coupling_diagnostics_driver_mod, only : multires_coupling_diagnostics_driver
  use field_collection_mod,                     only : field_collection_type
  use field_mod,                                only : field_type

  implicit none

  private
  public initialise, run, finalise

  type(model_clock_type), allocatable :: model_clock

  ! Dynamics and Physics mesh ids
  type(mesh_type), pointer :: dynamics_mesh    => null()
  type(mesh_type), pointer :: dynamics_2D_mesh => null()
  type(mesh_type), pointer :: physics_mesh     => null()
  type(mesh_type), pointer :: physics_2D_mesh  => null()
  type(mesh_type), pointer :: aerosol_mesh     => null()
  type(mesh_type), pointer :: aerosol_2D_mesh  => null()

  ! 2D Mesh names
  character(str_def) :: dynamics_2D_mesh_name
  character(str_def) :: physics_2D_mesh_name

  character(*), public, parameter :: xios_ctx = 'multires_coupling'

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !>
  !> @param [in,out] dynamics_mesh_modeldb The structure that holds model
  !>                                       state for the dynamics mesh
  !> @param [in,out] physics_mesh_modeldb  The structure that holds model
  !>                                       state for the physics mesh
  subroutine initialise( dynamics_mesh_modeldb, &
                         physics_mesh_modeldb,  &
                         program_name )

    implicit none

    type(modeldb_type), intent(inout) :: dynamics_mesh_modeldb
    type(modeldb_type), intent(inout) :: physics_mesh_modeldb
    character(*),    intent(in)    :: program_name

    !-------------------------------------------------------------------------
    ! Model init
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! Initialise Infrastructure such as meshes, FEM and runtime constants
    !-------------------------------------------------------------------------
    call log_event( 'Initialising Infrastructure...', LOG_LEVEL_ALWAYS )

    call initialise_infrastructure( program_name,            &
                                    model_clock,             &
                                    dynamics_mesh_modeldb%mpi )

    dynamics_2D_mesh_name = trim(dynamics_mesh_name)//'_2d'

    dynamics_mesh    => mesh_collection%get_mesh( dynamics_mesh_name )
    dynamics_2D_mesh => mesh_collection%get_mesh( dynamics_2D_mesh_name )

    physics_2D_mesh_name = trim(physics_mesh_name)//'_2d'
    physics_mesh    => mesh_collection%get_mesh( physics_mesh_name )
    physics_2D_mesh => mesh_collection%get_mesh( physics_2D_mesh_name )

    !-------------------------------------------------------------------------
    ! Create and initialise model data
    !-------------------------------------------------------------------------
    call log_event( 'Creating and Initialising Model Data...', LOG_LEVEL_ALWAYS )

    ! Instantiate the fields stored in model_data
    call create_model_data( dynamics_mesh_modeldb%model_data, &
                            dynamics_mesh,                    &
                            dynamics_2D_mesh,                 &
                            aerosol_mesh,                     &
                            aerosol_2D_mesh,                  &
                            model_clock )
    call create_model_data( physics_mesh_modeldb%model_data, &
                            physics_mesh,                    &
                            physics_2D_mesh,                 &
                            aerosol_mesh,                    &
                            aerosol_2D_mesh,                 &
                            model_clock )


    ! Initialise the fields stored in the model_data
    call initialise_model_data( dynamics_mesh_modeldb%model_data, &
                                model_clock,                      &
                                dynamics_mesh,                    &
                                dynamics_2D_mesh )
    call initialise_model_data( physics_mesh_modeldb%model_data,  &
                                model_clock,                      &
                                physics_mesh,                     &
                                physics_2D_mesh )

   ! Initial output
   ! We only want these once at the beginning of a run
   if ( model_clock%is_initialisation() .and. write_diag .and. &
        multires_coupling_mode /= multires_coupling_mode_test ) then
     call multires_coupling_diagnostics_driver( dynamics_mesh,            &
                                                dynamics_2D_mesh,         &
                                        dynamics_mesh_modeldb%model_data, &
                                                model_clock,              &
                                                nodal_output_on_w3 )
   end if

   if ( multires_coupling_mode /= multires_coupling_mode_test ) then
     ! Model configuration initialisation
     call initialise_model( dynamics_mesh,            &
                    dynamics_mesh_modeldb%model_data, &
                            physics_mesh,             &
                    physics_mesh_modeldb%model_data )
   end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>
  !> @brief Timesteps the model, calling the desired timestepping algorithm
  !> @param [in,out] dynamics_mesh_modeldb The structure that holds model
  !>                                       state for the dynamics mesh
  !> @param [in,out] physics_mesh_modeldb  The structure that holds model
  !>                                       state for the physics mesh
  subroutine run( dynamics_mesh_modeldb, physics_mesh_modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: dynamics_mesh_modeldb
    type(modeldb_type), intent(inout) :: physics_mesh_modeldb

    if ( multires_coupling_mode == multires_coupling_mode_test ) then
      ! Call coupling test algorithm
      call coupling_test_alg()
    else
      ! Do timestep
      do while (model_clock%tick())
        ! Update time-varying fields
        call update_variable_fields(                         &
          dynamics_mesh_modeldb%model_data%ancil_times_list, &
          model_clock,                                       &
          dynamics_mesh_modeldb%model_data%ancil_fields )
        ! Perform gungho timestep
        call gungho_step( dynamics_mesh,         &
                          dynamics_2D_mesh,      &
                          dynamics_mesh_modeldb, &
                          model_clock )
        ! Write out output file
        call log_event( "Writing depository output", LOG_LEVEL_INFO )

        if ( (mod(model_clock%get_step(), diagnostic_frequency) == 0) &
              .and. (write_diag) ) then
              call multires_coupling_diagnostics_driver( &
                dynamics_mesh,                    &
                dynamics_2D_mesh,                 &
                dynamics_mesh_modeldb%model_data, &
                model_clock,                      &
                nodal_output_on_w3 )
        end if
      end do
    end if

  end subroutine run

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Finalise after run
  !> @param [in,out] dynamics_mesh_modeldb The structure that holds model
  !>                                       state for the dynamics mesh
  !> @param [in,out] physics_mesh_modeldb  The structure that holds model
  !>                                       state for the physics mesh
  !>
  subroutine finalise( dynamics_mesh_modeldb, &
                       physics_mesh_modeldb,  &
                       program_name )

    implicit none

    type(modeldb_type), intent(inout) :: dynamics_mesh_modeldb
    type(modeldb_type), intent(inout) :: physics_mesh_modeldb

    character(*), intent(in) :: program_name

    !--------------------------------------------------------------------------
    ! Model finalise
    !--------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! Driver layer finalise
    !-------------------------------------------------------------------------

   ! Finalise XIOS context if we used it for diagnostic output or checkpointing
    if ( use_xios_io ) then
      call xios_context_finalize()
    end if

    if ( multires_coupling_mode /= multires_coupling_mode_test ) then
      call finalise_model( dynamics_mesh_modeldb%model_data, &
                           physics_mesh_modeldb%model_data,  &
                           program_name )
    end if

    ! Destroy the fields stored in model_data
    call finalise_model_data( dynamics_mesh_modeldb%model_data )
    call finalise_model_data( physics_mesh_modeldb%model_data )

    call finalise_infrastructure( program_name )

    call log_event( 'Miniapp completed', LOG_LEVEL_INFO )

  end subroutine finalise

end module multires_coupling_driver_mod
