!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> Drives the execution of the gravity_wave miniapp.
!>
module gravity_wave_driver_mod

  use constants_mod,                  only: i_def, i_native
  use gravity_wave_infrastructure_mod,only: initialise_infrastructure, &
                                            finalise_infrastructure
  use gravity_wave_grid_mod,          only: initialise_grid
  use gravity_wave_io_mod,            only: initialise_io, &
                                            finalise_io
  use create_gravity_wave_prognostics_mod, &
                                      only: create_gravity_wave_prognostics
  use runtime_constants_mod,          only: create_runtime_constants
  use gravity_wave_diagnostics_driver_mod, &
                                      only: gravity_wave_diagnostics_driver
  use field_mod,                      only: field_type
  use field_collection_mod,           only: field_collection_type
  use function_space_chain_mod,       only: function_space_chain_type
  use gravity_wave_mod,               only: program_name
  use gw_init_fields_alg_mod,         only: gw_init_fields_alg
  use init_gravity_wave_mod,          only: init_gravity_wave
  use step_gravity_wave_mod,          only: step_gravity_wave
  use final_gravity_wave_mod,         only: final_gravity_wave
  use log_mod,                        only: log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ALWAYS,   &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_TRACE
  use read_methods_mod,               only: read_checkpoint
  use write_methods_mod,              only: write_checkpoint
  use io_config_mod,                  only: write_diag,           &
                                            checkpoint_read,      &
                                            checkpoint_write,     &
                                            diagnostic_frequency, &
                                            use_xios_io,          &
                                            nodal_output_on_w3,   &
                                            subroutine_timers
  use files_config_mod,               only: checkpoint_stem_name
  use time_config_mod,                only: timestep_start, &
                                            timestep_end
  use timer_mod,                      only: init_timer, timer, output_timer
  use xios,                           only: xios_update_calendar

  implicit none

  private

  public initialise, run, finalise

  character(*), public, parameter   :: xios_ctx = 'gravity_wave'

  ! Depository of shared fields
  type( field_collection_type ), target :: depository

  ! Gravity wave model state
  type( field_collection_type ) :: prognostic_fields

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

  ! Function space chains
  type(function_space_chain_type) :: multigrid_function_space_chain

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

  implicit none

  character(*),      intent(in) :: filename
  integer(i_native), intent(in) :: model_communicator

  integer(i_def) :: ts_init

  ! Initialise aspects of the infrastructure
  call initialise_infrastructure( model_communicator, &
                                  filename,           &
                                  program_name )

  call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  if ( subroutine_timers ) then
    call init_timer()
    call timer(program_name)
  end if

  multigrid_function_space_chain = function_space_chain_type()

  ! Initialise aspects of the grid
  call initialise_grid(mesh_id, twod_mesh_id, chi, &
                       multigrid_function_space_chain)

  !Initialise aspects of output
  call initialise_io( model_communicator, &
                      mesh_id,            &
                      twod_mesh_id,       &
                      chi,                &
                      xios_ctx )

  ! Create runtime_constants object. This in turn creates various things
  ! needed by the timestepping algorithms such as mass matrix operators, mass
  ! matrix diagonal fields and the geopotential field
  call create_runtime_constants(mesh_id, twod_mesh_id, chi)

  ! Create the depository and prognostics field collections.
  ! Actually, here they will have the same contents (prognostics points to all
  ! the fields in the depository), but both are maintained to be consistent
  ! with more complex models
  depository=field_collection_type(name='depository')
  prognostic_fields=field_collection_type(name='prognostics')

  ! Create the prognostic fields
  call create_gravity_wave_prognostics(mesh_id, depository, prognostic_fields)

  ! Initialise prognostic fields
  if (checkpoint_read) then                 ! Recorded check point to start from
    call read_checkpoint(prognostic_fields, timestep_start-1)
  else                                      ! No check point to start from
    call gw_init_fields_alg(prognostic_fields)
  end if

  ! Initialise the gravity-wave model
  call init_gravity_wave( mesh_id, prognostic_fields )

  ! Output initial conditions
  ! We only want these once at the beginning of a run
  ts_init = max( (timestep_start - 1), 0 ) ! 0 or t previous.

  if (ts_init == 0) then

    if ( use_xios_io ) then

      ! Need to ensure calendar is initialised here as XIOS has no concept of timestep 0
      call xios_update_calendar(ts_init + 1)

    end if

    if ( write_diag ) then
      call gravity_wave_diagnostics_driver( mesh_id,           &
                                            prognostic_fields, &
                                            ts_init,           &
                                            nodal_output_on_w3)
    end if

  end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

  implicit none

  integer(i_def) :: timestep

  write(log_scratch_space,'(A,I0,A)') 'Running '//program_name//' ...'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  !--------------------------------------------------------------------------
  ! Model step
  !--------------------------------------------------------------------------
  do timestep = timestep_start,timestep_end

    ! Update XIOS calendar if we are using it for diagnostic output or checkpoint
    if ( use_xios_io ) then
      call log_event( program_name//': Updating XIOS timestep', LOG_LEVEL_INFO )
      call xios_update_calendar(timestep)
    end if

    call log_event( &
    "/****************************************************************************\ ", &
    LOG_LEVEL_TRACE )
    write( log_scratch_space, '(A,I0)' ) 'Start of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call step_gravity_wave(prognostic_fields)

    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( &
    '\****************************************************************************/ ', &
    LOG_LEVEL_INFO )
    if ( (mod(timestep, diagnostic_frequency) == 0) .and. (write_diag) ) then

      call log_event("Gravity Wave: writing diagnostic output", LOG_LEVEL_INFO)

      call gravity_wave_diagnostics_driver( mesh_id,           &
                                            prognostic_fields, &
                                            timestep,          &
                                            nodal_output_on_w3)
    end if

  end do

  end subroutine run


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

  implicit none

  !--------------------------------------------------------------------------
  ! Model finalise
  !--------------------------------------------------------------------------
  call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  call final_gravity_wave(prognostic_fields, program_name)

  ! Write checkpoint/restart files if required
  if( checkpoint_write ) then
    call write_checkpoint(prognostic_fields,timestep_end)
  end if

  if ( subroutine_timers ) then
    call timer(program_name)
    call output_timer()
  end if

  !--------------------------------------------------------------------------
  ! Driver layer finalise
  !--------------------------------------------------------------------------

  call finalise_io()

  call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

  call finalise_infrastructure()

  end subroutine finalise

end module gravity_wave_driver_mod
