!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @mainpage solver_miniapp
!> The solver API is abstract so cannot be tested unless with a particular implementation.
!! The mini-app makes an example linear operator, preconditioner and krylov subspace solver.
!! This can then be executed to test the solver api against different compilers.

program solver_miniapp

  use constants_mod,                    only : i_def, PRECISION_REAL
  use convert_to_upper_mod,             only : convert_to_upper
  use cli_mod,                          only : get_initial_filename
  use driver_mesh_mod,                  only : init_mesh
  use driver_fem_mod,                   only : init_fem
  use driver_log_mod,                   only : init_logger, final_logger
  use halo_comms_mod,                   only : initialise_halo_comms, &
                                               finalise_halo_comms
  use init_solver_miniapp_mod,          only : init_solver_miniapp
  use mpi_mod,                          only : initialise_comm, store_comm, &
                                               finalise_comm,               &
                                               get_comm_size, get_comm_rank
  use field_mod,                        only : field_type
  use field_vector_mod,                 only : field_vector_type
  use solver_miniapp_alg_mod,           only : solver_miniapp_alg
  use configuration_mod,                only : final_configuration
  use solver_miniapp_mod,               only : load_configuration, program_name
  use log_mod,                          only : log_event,            &
                                               log_scratch_space,    &
                                               LOG_LEVEL_ALWAYS,     &
                                               LOG_LEVEL_INFO
  use mesh_mod,                         only : mesh_type
  use checksum_alg_mod,                 only : checksum_alg

  implicit none

  character(:), allocatable :: filename

  integer(i_def) :: total_ranks, local_rank
  integer(i_def) :: comm = -999

  ! prognostic fields
  type( field_type ), target, dimension(3) :: chi
  type( field_type ), target :: panel_id
  type( field_type ) :: field_1, field_2
  type( field_vector_type) :: fv_1

  type(mesh_type), pointer :: mesh => null()


  !-----------------------------------------------------------------------------
  ! Driver layer init
  !-----------------------------------------------------------------------------

  ! Initialise MPI communicatios and get a valid communicator
  call initialise_comm(comm)

  ! Initialise halo functionality
  call initialise_halo_comms( comm )

  ! Save the commmunicator for later use
  call store_comm(comm)

  total_ranks = get_comm_size()
  local_rank  = get_comm_rank()

  call get_initial_filename( filename )
  call load_configuration( filename )
  deallocate( filename )

  call init_logger(total_ranks, local_rank, program_name)

  write(log_scratch_space,'(A)')                        &
      'Application built with '//trim(PRECISION_REAL)// &
      '-bit real numbers'
  call log_event( log_scratch_space, LOG_LEVEL_ALWAYS )

  !-----------------------------------------------------------------------------
  ! model init
  !-----------------------------------------------------------------------------
  call log_event( 'Initialising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  call init_mesh( local_rank, total_ranks, mesh )

  call init_fem( mesh, chi, panel_id )

  ! Create and initialise prognostic fields
  call init_solver_miniapp( mesh, chi, panel_id, fv_1 )

  ! Call an algorithm
  call solver_miniapp_alg( fv_1, chi, panel_id )

  ! Write out output file
  call log_event(program_name//": writing diagnostic output", LOG_LEVEL_INFO)

  ! pull the fields from the vector
  call fv_1%export_field( field_1, 1 )
  call fv_1%export_field( field_2, 2 )

  !-----------------------------------------------------------------------------
  ! model finalise
  !-----------------------------------------------------------------------------
  call log_event( 'Finalising '//program_name//' ...', LOG_LEVEL_ALWAYS )

  ! Write checksums to file
  call checksum_alg(program_name, field_1, 'solver_field_1',field_2, 'solver_field_2')

  call log_event( program_name//': Miniapp completed', LOG_LEVEL_INFO )

  !-----------------------------------------------------------------------------
  ! Driver layer finalise
  !-----------------------------------------------------------------------------

  ! Finalise namelist configurations
  call final_configuration()

  ! Finalise halo functionality
  call finalise_halo_comms()

  ! Finalise MPI communications
  call finalise_comm()

  ! Finalise the logging system
  call final_logger(program_name)

end program solver_miniapp
