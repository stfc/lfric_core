!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the io_dev miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module io_dev_driver_mod

use configuration_mod,          only: final_configuration
use constants_mod,              only: i_def, i_native, PRECISION_REAL
use convert_to_upper_mod,       only: convert_to_upper
use create_fem_mod,             only: init_fem
use create_mesh_mod,            only: init_mesh
use derived_config_mod,         only: set_derived_config
use diagnostics_io_mod,         only: write_scalar_diagnostic, &
                                      write_vector_diagnostic
use field_mod,                  only: field_type
use finite_element_config_mod,  only: element_order
use global_mesh_collection_mod, only: global_mesh_collection, &
                                      global_mesh_collection_type
use init_io_dev_mod,            only: init_io_dev
use io_dev_mod,                 only: load_configuration, program_name
use io_mod,                     only: initialise_xios
use log_mod,                    only: log_event,          &
                                      log_set_level,      &
                                      log_scratch_space,  &
                                      initialise_logging, &
                                      finalise_logging,   &
                                      LOG_LEVEL_ALWAYS,   &
                                      LOG_LEVEL_ERROR,    &
                                      LOG_LEVEL_WARNING,  &
                                      LOG_LEVEL_INFO,     &
                                      LOG_LEVEL_DEBUG,    &
                                      LOG_LEVEL_TRACE
use mpi_mod,                    only: store_comm,    &
                                      get_comm_size, &
                                      get_comm_rank
use timestepping_config_mod,    only: dt
use xios,                       only: xios_context_finalize, &
                                      xios_update_calendar
use yaxt,                       only: xt_initialize, xt_finalize

implicit none

private

public initialise, run, finalise

character(*), public, parameter :: xios_ctx = 'io_dev'

! Examples of fields that are output on different I/O domains
type( field_type ) :: density ! on W3
type( field_type ) :: theta   ! on Wtheta
type( field_type ) :: wind    ! on W2 and output on W2H and W2V


! Coordinate field
type(field_type), target, dimension(3) :: chi

integer(i_def) :: mesh_id
integer(i_def) :: twod_mesh_id

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>@param[in] filename Name of the file containing the desired configuration
  subroutine initialise( filename, model_communicator )

  use logging_config_mod, only: run_log_level,          &
                                key_from_run_log_level, &
                                RUN_LOG_LEVEL_ERROR,    &
                                RUN_LOG_LEVEL_INFO,     &
                                RUN_LOG_LEVEL_DEBUG,    &
                                RUN_LOG_LEVEL_TRACE,    &
                                RUN_LOG_LEVEL_WARNING

  implicit none

  character(*),      intent(in) :: filename
  integer(i_native), intent(in) :: model_communicator

  integer(i_def) :: total_ranks, local_rank

  integer(i_def) :: dtime
  integer(i_native) :: log_level

  ! Save lfric's part of the split communicator for later use
  call store_comm( model_communicator )

  ! Initialise YAXT
  call xt_initialize( model_communicator )

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call initialise_logging(local_rank, total_ranks, program_name)

  call load_configuration( filename )

  select case (run_log_level)
  case( RUN_LOG_LEVEL_ERROR )
    log_level = LOG_LEVEL_ERROR
  case( RUN_LOG_LEVEL_WARNING )
    log_level = LOG_LEVEL_WARNING
  case( RUN_LOG_LEVEL_INFO )
    log_level = LOG_LEVEL_INFO
  case( RUN_LOG_LEVEL_DEBUG )
    log_level = LOG_LEVEL_DEBUG
  case( RUN_LOG_LEVEL_TRACE )
    log_level = LOG_LEVEL_TRACE
  end select

  call log_set_level( log_level )

  write(log_scratch_space,'(A)')                              &
      'Runtime message logging severity set to log level: '// &
      convert_to_upper(key_from_run_log_level(run_log_level))
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  write(log_scratch_space,'(A)')                        &
      'Application built with '//trim(PRECISION_REAL)// &
      '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  call set_derived_config( .true. )

  !----------------------------------------------------------------------------
  ! Mesh init
  !----------------------------------------------------------------------------
  call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  allocate( global_mesh_collection, &
            source = global_mesh_collection_type() )

  ! Create the mesh
  call init_mesh(local_rank, total_ranks, mesh_id, twod_mesh_id)

  !----------------------------------------------------------------------------
  ! FEM init
  !----------------------------------------------------------------------------
  ! Create FEM specifics (function spaces and chi field)
  call init_fem(mesh_id, chi)


  !-----------------------------------------------------------------------------
  ! IO init
  !-----------------------------------------------------------------------------

  ! Set up XIOS domain and context
  dtime = int(dt)

  call initialise_xios( xios_ctx,     &
                        model_communicator, &
                        dtime,        &
                        mesh_id,      &
                        twod_mesh_id, &
                        chi )

  !-----------------------------------------------------------------------------
  ! Model init
  !-----------------------------------------------------------------------------

  ! Create and initialise fields
  call init_io_dev(mesh_id, twod_mesh_id, chi, density, theta, wind)


  ! Full global meshes no longer required, so reclaim
  ! the memory from global_mesh_collection
  write(log_scratch_space,'(A)') &
      "Purging global mesh collection."
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  deallocate(global_mesh_collection)

  end subroutine initialise


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Timesteps the model, calling the desired timestepping algorithm based
  !>upon the configuration
  subroutine run()

  implicit none

  write(log_scratch_space,'(A)') 'Running '//program_name//' ...'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  !-----------------------------------------------------------------------------
  ! Model step
  !-----------------------------------------------------------------------------

  ! Update XIOS calendar
  call log_event( program_name//': Updating XIOS timestep', LOG_LEVEL_INFO )
  call xios_update_calendar(1)

  ! Write out the fields

  call log_event( program_name//': Writing XIOS output', LOG_LEVEL_INFO)
  call write_scalar_diagnostic('density', density, 1, mesh_id, .false.)
  call write_scalar_diagnostic('theta', theta, 1, mesh_id, .false.)
  call write_vector_diagnostic('wind', wind, 1, mesh_id, .false.)

  end subroutine run

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !>@brief Tidies up after a run.

  subroutine finalise()

  implicit none

  call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  ! Finalise XIOS context
  call xios_context_finalize()

  ! Finalise namelist configurations
  call final_configuration()

  ! Finalise YAXT
  call xt_finalize()

  call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

  ! Finalise the logging system
  call finalise_logging()

  end subroutine finalise

end module io_dev_driver_mod
