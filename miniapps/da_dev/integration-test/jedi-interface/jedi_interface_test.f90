!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the da dev, jedi interface
!!         integrationf tests.
!>@details Sets up and runs the integration tests specified in
!!         jedi_interface_test.py.
program jedi_interface_test

  use configuration_mod,              only : final_configuration, &
                                             read_configuration
  use constants_mod,                  only : i_def, r_def, l_def
  use halo_comms_mod,                 only : initialise_halo_comms, &
                                             finalise_halo_comms
  use test_jedi_interface_driver_mod, only : test_jedi_interface_init,          &
                                             test_jedi_interface_final,         &
                                             run_init_lfric_calendar_start,     &
                                             run_init_lfric_calendar_start_err, &
                                             run_init_string_err,               &
                                             run_init_from_jedi_datetime_err,   &
                                             run_YYYYMMDD_to_JDN,               &
                                             run_JDN_to_YYYYMMDD_invalid,       &
                                             run_hhmmss_to_seconds,             &
                                             run_seconds_to_hhmmss_large,       &
                                             run_seconds_to_hhmmss_neg
  use log_mod,                        only : log_event,          &
                                             initialise_logging, &
                                             finalise_logging,   &
                                             LOG_LEVEL_ERROR,    &
                                             LOG_LEVEL_INFO
  use mpi_mod,                        only : global_mpi, &
                                             create_comm, destroy_comm

  implicit none

  ! MPI communicator
  integer(i_def) :: comm

  ! Number of processes and local rank
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer(i_def)   :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag

  ! Flags which determine the tests that will be carried out
  logical(l_def) :: do_test_init_lfric_calendar_start = .false.
  logical(l_def) :: do_test_init_lfric_calendar_start_err = .false.
  logical(l_def) :: do_test_init_string_err = .false.
  logical(l_def) :: do_test_init_from_jedi_datetime_err = .false.
  logical(l_def) :: do_test_YYYYMMDD_to_JDN = .false.
  logical(l_def) :: do_test_JDN_to_YYYYMMDD_invalid = .false.
  logical(l_def) :: do_test_hhmmss_to_seconds = .false.
  logical(l_def) :: do_test_seconds_to_hhmmss_large = .false.
  logical(l_def) :: do_test_seconds_to_hhmmss_neg = .false.

  ! Usage message to print
  character(len=256) :: usage_message

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

  call initialise_logging( comm, 'jedi-interface_test' )

  call log_event( 'jedi interface testing running ...', LOG_LEVEL_INFO )

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
    write(usage_message,*) "Usage: ",trim(program_name), &
      " <namelist filename> "            // &
      " test_XXX with XXX in { "         // &
      " init_lfric_calendar_start, "     // &
      " init_lfric_calendar_start_err, " // &
      " init_string_err, "               // &
      " init_from_jedi_datetime_err, "   // &
      " YYYYMMDD_to_JDN, "               // &
      " JDN_to_YYYYMMDD_invalid, "       // &
      " hhmmss_to_seconds, "             // &
      " seconds_to_hhmmss_large, "       // &
      " seconds_to_hhmmss_neg, "         // &
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
  case ("test_init_lfric_calendar_start")
    do_test_init_lfric_calendar_start = .true.
  case ("test_init_lfric_calendar_start_err")
    do_test_init_lfric_calendar_start_err = .true.
  case ("test_init_string_err")
    do_test_init_string_err = .true.
  case ("test_init_from_jedi_datetime_err")
    do_test_init_from_jedi_datetime_err = .true.
  case ("test_YYYYMMDD_to_JDN")
    do_test_YYYYMMDD_to_JDN = .true.
  case ("test_JDN_to_YYYYMMDD_invalid")
    do_test_JDN_to_YYYYMMDD_invalid = .true.
  case ("test_hhmmss_to_seconds")
    do_test_hhmmss_to_seconds = .true.
  case ("test_seconds_to_hhmmss_large")
    do_test_seconds_to_hhmmss_large = .true.
  case ("test_seconds_to_hhmmss_neg")
    do_test_seconds_to_hhmmss_neg = .true.
  case default
    call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  ! must be before load_configuration call
  if ( do_test_init_lfric_calendar_start_err ) then
    call run_init_lfric_calendar_start_err()
  end if

  ! Setup configuration, and initialise tests
  call read_configuration( filename )
  call test_jedi_interface_init()

  if ( do_test_init_lfric_calendar_start ) then
    call run_init_lfric_calendar_start()
  end if
  if ( do_test_init_string_err ) then
    call run_init_string_err()
  end if
  if ( do_test_init_from_jedi_datetime_err ) then
    call run_init_from_jedi_datetime_err()
  end if
  if ( do_test_YYYYMMDD_to_JDN ) then
    call run_YYYYMMDD_to_JDN()
  end if
  if ( do_test_JDN_to_YYYYMMDD_invalid ) then
    call run_JDN_to_YYYYMMDD_invalid()
  end if
  if ( do_test_hhmmss_to_seconds ) then
    call run_hhmmss_to_seconds()
  end if
  if ( do_test_seconds_to_hhmmss_large ) then
    call run_seconds_to_hhmmss_large()
  end if
  if ( do_test_seconds_to_hhmmss_neg ) then
    call run_seconds_to_hhmmss_neg()
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down

  call test_jedi_interface_final()

  if (allocated(program_name)) deallocate(program_name)
  if (allocated(filename))     deallocate(filename)
  if (allocated(test_flag))    deallocate(test_flag)

  call log_event( 'jedi-interface functional testing completed ...', LOG_LEVEL_INFO )

  ! Finalise halo functionality
  call finalise_halo_comms()

  ! Finalise the logging system
  call global_mpi%finalise()

  call final_configuration

  call destroy_comm()

end program jedi_interface_test
