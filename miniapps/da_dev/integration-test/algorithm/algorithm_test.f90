!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the da dev algorithm
!!         integration tests.
!>@details Sets up and runs the integration tests specified in
!!         algorithms_test.py. Currently only
!!         da_dev_increment_alg_mod.x90 is tested, this algorithm
!!         takes a real field and adds one to its value.
program algorithm_test

  use configuration_mod,             only : final_configuration, &
                                            read_configuration
  use constants_mod,                 only : i_def, r_def
  use test_algorithm_mod,            only : test_algorithm_finalise,   &
                                            test_algorithm_initialise, &
                                            test_da_dev_increment_alg
  use driver_mesh_mod,               only : init_mesh, final_mesh
  use halo_comms_mod,                only : initialise_halo_comms, &
                                            finalise_halo_comms
  use log_mod,                       only : log_event,          &
                                            initialise_logging, &
                                            finalise_logging,   &
                                            LOG_LEVEL_ERROR,    &
                                            LOG_LEVEL_INFO
  use mesh_mod,                      only : mesh_type
  use mpi_mod,                       only : global_mpi, &
                                            create_comm, destroy_comm

  implicit none

  ! MPI communicator
  integer(i_def) :: comm

  ! Number of processes and local rank
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag

  ! Flags which determine the tests that will be carried out
  logical :: do_test_da_dev_increment_alg_mod = .false.

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Setup for all tests
  ! Usage message to print
  character(len=256)       :: usage_message

  type(mesh_type), pointer :: mesh => null()

  real(r_def), parameter   :: tolerance = 1.0e-3_r_def

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Communicators and Logging Setup
  ! Initialise MPI communicatios and get a valid communicator
  call create_comm(comm)

  ! Save the communicator for later use
  call global_mpi%initialise(comm)

  ! Initialise halo functionality
  call initialise_halo_comms(comm)

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()

  call initialise_logging( comm, 'da_dev_test_alg' )

  call log_event( 'da dev alg testing running ...', LOG_LEVEL_INFO )

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
    write(usage_message,*) "Usage: ",trim(program_name), &
      " <namelist filename> "       // &
      " test_XXX with XXX in { "    // &
      " da_dev_increment_alg_mod, " // &
      " } "
    call log_event( trim(usage_message), LOG_LEVEL_ERROR )
  end if

  call get_command_argument( 1, dummy, length, status )
  allocate( character(length) :: filename )
  call get_command_argument( 1, filename, length, status )

  call get_command_argument( 2, dummy, length, status )
  allocate(character(length)::test_flag)
  call get_command_argument( 2, test_flag, length, status )

  ! Choose test case depending on flag provided in the first command
  ! line argument
  select case (trim(test_flag))
  case ("test_da_dev_increment_alg_mod")
    do_test_da_dev_increment_alg_mod = .true.
  case default
    call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  ! Setup configuration, mesh, and fem
  call read_configuration( filename )
  call init_mesh( local_rank, total_ranks, mesh )
  call test_algorithm_initialise(mesh) ! fem

  if ( do_test_da_dev_increment_alg_mod ) then
    call test_da_dev_increment_alg(tolerance)
  endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down

  call test_algorithm_finalise()
  call final_mesh()
  nullify(mesh)

  if (allocated(program_name)) deallocate(program_name)
  if (allocated(filename))     deallocate(filename)
  if (allocated(test_flag))    deallocate(test_flag)

  call log_event( 'da dev alg functional testing completed ...', LOG_LEVEL_INFO )

  ! Finalise halo functionality
  call finalise_halo_comms()

  ! Finalise MPI communications
  call global_mpi%finalise()
  call destroy_comm()

  ! Finalise the logging system
  call finalise_logging()

  call final_configuration()

end program algorithm_test
