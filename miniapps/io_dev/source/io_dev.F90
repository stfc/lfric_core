!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page io_dev io_dev Miniapp
!> Test program for the XIOS IO implementation.

!> @brief Main program used to test XIOS setup and output of a field.

program io_dev

  use cli_mod,                only: get_initial_filename
  use constants_mod,          only: precision_real
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use driver_time_mod,        only: init_time, final_time
  use driver_timer_mod,       only: init_timers, final_timers
  use io_dev_mod,             only: io_dev_required_namelists
  use io_dev_driver_mod,      only: initialise, step, finalise
  use io_dev_modeldb_mod,     only: modeldb_type
  use log_mod,                only: log_event,       &
                                    log_level_trace, &
                                    log_scratch_space
  use model_clock_mod,        only: model_clock_type
  use mpi_mod,                only: global_mpi

  use namelist_collection_mod, only: namelist_collection_type

  implicit none

  character(*), parameter :: program_name = "io_dev"

  type (modeldb_type)             :: modeldb
  character(:),       allocatable :: filename

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( program_name, table_len=10 )

  call modeldb%fields%add_empty_field_collection("depository", table_len=1)
  call modeldb%fields%add_empty_field_collection("dump_fields", table_len=1)
  call modeldb%fields%add_empty_field_collection("alg_fields", table_len=1)

  call modeldb%io_contexts%initialise(program_name, 100)

  write( log_scratch_space,'(A)' )                         &
      'Application built with ' // trim(precision_real) // &
      '-bit real numbers.'
  call log_event( log_scratch_space, log_level_trace )

  call init_comm( "io_dev", global_mpi )
  call get_initial_filename( filename )
  call init_config( filename, io_dev_required_namelists, &
                    modeldb%configuration )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_collections()
  call init_time( modeldb%clock, modeldb%calendar )
  deallocate( filename )

  call log_event( 'Initialising '//program_name//' ...', log_level_trace )
  call initialise( program_name, modeldb )

  write(log_scratch_space,'("Running ", A, " ...")') program_name
  call log_event( log_scratch_space, log_level_trace )
  do while( modeldb%clock%tick() )
    call step( modeldb, program_name )
  end do

  call log_event( 'Finalising '//program_name//' ...', log_level_trace )
  call finalise( modeldb )

  call final_time( modeldb%clock, modeldb%calendar )
  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( global_mpi )

end program io_dev
