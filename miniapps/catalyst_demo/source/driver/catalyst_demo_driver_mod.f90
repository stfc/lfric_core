!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the catalyst_demo miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module catalyst_demo_driver_mod

  use checksum_alg_mod,               only: checksum_alg
  use constants_mod,                  only: i_def, i_native
  use derived_config_mod,             only: set_derived_config
  use yaxt,                           only: xt_initialize, xt_finalize
  use field_mod,                      only: field_type
  use finite_element_config_mod,      only: element_order
  use function_space_chain_mod,       only: function_space_chain_type
  use global_mesh_collection_mod,     only: global_mesh_collection, &
                                            global_mesh_collection_type
  use catalyst_demo_mod,              only: load_configuration
  use gw_alg_mod,                     only: gravity_wave_alg_init, &
                                            gravity_wave_alg_step, &
                                            gravity_wave_alg_final
  use create_fem_mod,                 only: init_fem
  use create_mesh_mod,                only: init_mesh
  use init_catalyst_demo_mod,         only: init_catalyst_demo
  use io_mod,                         only: initialise_xios,   &
                                            ts_fname
  use diagnostics_io_mod,             only: write_scalar_diagnostic, &
                                            write_vector_diagnostic
  use log_mod,                        only: log_event,          &
                                            log_set_level,      &
                                            log_scratch_space,  &
                                            initialise_logging, &
                                            finalise_logging,   &
                                            LOG_LEVEL_ERROR,    &
                                            LOG_LEVEL_INFO,     &
                                            LOG_LEVEL_DEBUG,    &
                                            LOG_LEVEL_TRACE,    &
                                            log_scratch_space
  use operator_mod,                   only: operator_type
  use io_config_mod,                  only: write_diag,           &
                                            diagnostic_frequency, &
                                            use_xios_io,          &
                                            nodal_output_on_w3,   &
                                            checkpoint_write,     &
                                            checkpoint_stem_name, &
                                            subroutine_timers
  use time_config_mod,                only: timestep_start, &
                                            timestep_end
  use timer_mod,                      only: init_timer, timer, output_timer
  use timestepping_config_mod,        only: dt
  use mpi_mod,                        only: store_comm,    &
                                            get_comm_size, &
                                            get_comm_rank
  use visualisation_config_mod,       only: write_catalyst_output,   &
                                            visualisation_frequency, &
                                            visualisation_stem_name, &
                                            use_python_vis_pipeline, &
                                            python_vis_pipeline_name
  use visualisation_mod,              only: vis_field_list_type, &
                                            catalyst_initialize, &
                                            catalyst_coprocess,  &
                                            catalyst_finalize
  use xios,                           only: xios_ctx,              &
                                            xios_id,               &
                                            xios_context_finalize, &
                                            xios_domain_init,      &
                                            xios_finalize,         &
                                            xios_initialize,       &
                                            xios_update_calendar
                                            initialise_xios,       &
                                            xios_finalize,         &
                                            xios_initialize,       &
                                            xios_update_calendar

  implicit none

  private

  public initialise, run, finalise

  character(*), public, parameter   :: program_name = 'catalyst_demo'
  character(*), public, parameter   :: xios_ctx = 'catalyst_demo'

  type(restart_type) :: restart

  ! Prognostic fields
  type( field_type ) :: wind
  type( field_type ) :: buoyancy
  type( field_type ) :: pressure

  ! Coordinate field
  type(field_type), target, dimension(3) :: chi

  integer(i_def) :: mesh_id
  integer(i_def) :: twod_mesh_id

  ! Function space chains
  type(function_space_chain_type) :: multigrid_function_space_chain

  ! Lists of fields for visualisation
  type(vis_field_list_type) :: vis_fields

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  subroutine initialise( filename, model_communicator )

  implicit none

  character(*),      intent(in) :: filename
  integer(i_native), intent(in) :: model_communicator

  integer(i_def) :: total_ranks, local_rank

  integer(i_def) :: ts_init
  integer(i_def) :: dtime

  !Store the MPI communicator for later use
  call store_comm(model_communicator)

  ! Initialise YAXT
  call xt_initialize(model_communicator)

  ! Get the rank information
  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, program_name)

  call log_event( program_name//': Running miniapp ...', LOG_LEVEL_INFO )

  call load_configuration( filename )
  call set_derived_config( .false. )

  restart = restart_type( restart_filename, local_rank, total_ranks )

  !----------------------------------------------------------------------------
  ! Mesh init
  !----------------------------------------------------------------------------
  if ( subroutine_timers ) then
    call init_timer()
    call timer(program_name)
  end if

  allocate( global_mesh_collection, &
       source = global_mesh_collection_type() )

  ! Create the mesh
  call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id )

  !----------------------------------------------------------------------------
  ! FEM init
  !----------------------------------------------------------------------------
  ! Create FEM specifics (function spaces and chi field)
  call init_fem(mesh_id, chi)

  !----------------------------------------------------------------------------
  ! IO init
  !----------------------------------------------------------------------------

  ! If using XIOS for diagnostic output or checkpointing, then set up XIOS
  ! domain and context
  if ( use_xios_io ) then

    dtime = int(dt)

    call initialise_xios( xios_ctx,     &
                           model_communicator, &
                          dtime,        &
                          mesh_id,      &
                          twod_mesh_id, &
                          chi)

    ! Make sure XIOS calendar is set to timestep 1 as it starts there
    ! not timestep 0.
    call xios_update_calendar(1)

  end if

  !----------------------------------------------------------------------------
  ! Visualisation init
  !----------------------------------------------------------------------------

  if ( write_catalyst_output ) then
    call catalyst_initialize(visualisation_frequency,       &
                             trim(visualisation_stem_name), &
                             model_communicator,            &
                             use_python_vis_pipeline,       &
                             trim(python_vis_pipeline_name))
  end if

  !----------------------------------------------------------------------------
  ! Model init
  !----------------------------------------------------------------------------
  multigrid_function_space_chain = function_space_chain_type()

  ! Create function space collection and initialise prognostic fields
  call init_catalyst_demo( mesh_id, twod_mesh_id, chi, multigrid_function_space_chain, &
                           wind, pressure, buoyancy )

  ! Full global meshes no longer required, so reclaim
  ! the memory from global_mesh_collection
  write(log_scratch_space,'(A)') &
      "Purging global mesh collection."
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  deallocate(global_mesh_collection)

  ! Output initial conditions
  ! We only want these once at the beginning of a run
  ts_init = max( (restart%ts_start() - 1), 0 ) ! 0 or t previous.

  if (ts_init == 0) then

    if ( use_xios_io ) then

      ! Need to ensure calendar is initialised here as XIOS has no concept of timestep 0
      call xios_update_calendar(ts_init + 1)

    end if

    if ( write_diag ) then
      ! Calculation and output of diagnostics
      call write_vector_diagnostic('wind', wind, ts_init, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('pressure', pressure, ts_init, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('buoyancy', buoyancy, ts_init, mesh_id, nodal_output_on_w3)
    end if

    ! Catalyst visualisation
    if (write_catalyst_output) then
      vis_fields = vis_field_list_type()
      call vis_fields%add_vis_field(wind, 'wind')
      call vis_fields%add_vis_field(pressure, 'pressure')
      call vis_fields%add_vis_field(buoyancy, 'buoyancy')
      call catalyst_coprocess(ts_init, ts_init*dt, vis_fields, mesh_id)
    end if

  end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time steps.
  !>
  subroutine run()

  implicit none

  integer(i_def) :: timestep

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
    if (timestep == timestep_start) then
      call gravity_wave_alg_init(mesh_id, wind, pressure, buoyancy)
    end if

    call gravity_wave_alg_step(wind, pressure, buoyancy)
    write( log_scratch_space, '(A,I0)' ) 'End of timestep ', timestep
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    call log_event( &
    '\****************************************************************************/ ', &
    LOG_LEVEL_INFO )
    if ( (mod(timestep, diagnostic_frequency) == 0) .and. (write_diag) ) then

      call log_event("Catalyst demo: writing diagnostic output", LOG_LEVEL_INFO)

      call write_vector_diagnostic('wind', wind, timestep, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('pressure', pressure, timestep, mesh_id, nodal_output_on_w3)
      call write_scalar_diagnostic('buoyancy', buoyancy, timestep, mesh_id, nodal_output_on_w3)

    end if

    ! Catalyst uses its own mechanism to determine if output should be written
    if (write_catalyst_output) then
      call catalyst_coprocess(timestep, timestep*dt, vis_fields, mesh_id)
    end if

  end do

  ! Write checkpoint/restart files if required
  if( checkpoint_write ) then
    write(log_scratch_space,'(A,I6)') &
        "Checkpointing pressure at timestep ", timestep_end
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call pressure%write_checkpoint("checkpoint_pressure", trim(ts_fname(checkpoint_stem_name,&
                                  "", "pressure", timestep_end,"")))

    write(log_scratch_space,'(A,I6)') &
         "Checkpointing wind at timestep ", timestep_end
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call wind%write_checkpoint("checkpoint_wind", trim(ts_fname(checkpoint_stem_name,&
                                  "", "wind", timestep_end,"")))

    write(log_scratch_space,'(A,I6)') &
         "Checkpointing buoyancy at timestep ", timestep_end
    call log_event(log_scratch_space,LOG_LEVEL_INFO)
    call buoyancy%write_checkpoint("checkpoint_buoyancy", trim(ts_fname(checkpoint_stem_name,&
                                  "", "buoyancy", timestep_end,"")))
  end if

  end subroutine run


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !>
  subroutine finalise()

  implicit none

  integer(i_def) :: rc
  integer(i_def) :: ierr

  !----------------------------------------------------------------------------
  ! Model finalise
  !----------------------------------------------------------------------------

  ! Write checksums to file
  call checksum_alg( program_name, wind, 'wind', buoyancy, 'buoyancy', pressure, 'pressure')

  call log_event( program_name//': Miniapp run complete', LOG_LEVEL_INFO )
  if ( subroutine_timers ) then
    call timer(program_name)
    call output_timer()
  end if

  call gravity_wave_alg_final()

  !----------------------------------------------------------------------------
  ! Driver layer finalise
  !----------------------------------------------------------------------------

  if ( write_catalyst_output ) then
    call catalyst_finalize()
  end if

  ! Finalise XIOS context if we used it for diagnostic output or checkpointing
  if ( use_xios_io ) then
    call xios_context_finalize()
  end if

  ! Finalise YAXT
  call xt_finalize()

  ! Finalise the logging system
  call finalise_logging()

end subroutine finalise

end module catalyst_demo_driver_mod
