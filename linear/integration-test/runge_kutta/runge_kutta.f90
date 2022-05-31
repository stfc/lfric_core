!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the tangent linear tests.
!>@details The default is to run all available tests - which
!!         test whether the linear code is tangent linear to the
!!         corresponding nonlinear code.
program runge_kutta

  use mod_wait,           only : init_wait
  use mpi_mod,            only : initialise_comm, &
                                 finalise_comm
  use xios,               only : xios_initialize
  use log_mod,            only : log_event,       &
                                 LOG_LEVEL_ERROR, &
                                 LOG_LEVEL_INFO
  use tl_test_driver_mod, only : initialise,                  &
                                 finalise,                    &
                                 run_timesteps,               &
                                 run_kinetic_energy_gradient, &
                                 run_advect_density_field,    &
                                 run_advect_theta_field,      &
                                 run_vorticity_advection,     &
                                 run_project_eos_pressure,    &
                                 run_hydrostatic,             &
                                 run_pressure_gradient_bd,    &
                                 run_rk_alg
  implicit none

  character(*), parameter :: xios_id = "linear"

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999
  integer                   :: model_communicator = -999

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag

  ! Flags which determine the tests that will be carried out
  logical :: do_test_timesteps = .false.
  logical :: do_test_kinetic_energy_gradient = .false.
  logical :: do_test_advect_density_field = .false.
  logical :: do_test_advect_theta_field = .false.
  logical :: do_test_project_eos_pressure = .false.
  logical :: do_test_vorticity_advection = .false.
  logical :: do_test_pressure_gradient_bd = .false.
  logical :: do_test_hydrostatic = .false.
  logical :: do_test_rk_alg = .false.

  ! Usage message to print
  character(len=256) :: usage_message

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  ! Initialise XIOS and get back the split mpi communicator. This requires
  ! an iodef.xml file to be available (even if use_xios_io is .false.).
  call init_wait()
  call xios_initialize(xios_id, return_comm = model_communicator)

  call log_event( 'TL testing running ...', LOG_LEVEL_INFO )

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
     write(usage_message,*) "Usage: ",trim(program_name), &
          " <namelist filename> "      // &
          " test_XXX with XXX in { "   // &
          " timesteps, "               // &
          " kinetic_energy_gradient, " // &
          " advect_density_field, "    // &
          " advect_theta_field, "      // &
          " vorticity_advection, "     // &
          " pressure_gradient_bd, "    // &
          " project_eos_pressure, "    // &
          " hydrostatic, "             // &
          " rk_alg, "                  // &
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
  case ("test_timesteps")
     do_test_timesteps = .true.
  case ("test_kinetic_energy_gradient")
     do_test_kinetic_energy_gradient = .true.
  case ("test_advect_density_field")
     do_test_advect_density_field = .true.
  case ("test_advect_theta_field")
     do_test_advect_theta_field = .true.
  case ("test_project_eos_pressure")
     do_test_project_eos_pressure = .true.
  case ("test_vorticity_advection")
     do_test_vorticity_advection = .true.
  case ("test_pressure_gradient_bd")
     do_test_pressure_gradient_bd = .true.
  case ("test_hydrostatic")
     do_test_hydrostatic = .true.
  case ("test_rk_alg")
     do_test_rk_alg = .true.
  case default
     call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  call initialise( filename, model_communicator )
  deallocate( filename )

  if (do_test_timesteps) then
    call run_timesteps()
  endif
  if (do_test_kinetic_energy_gradient) then
    call run_kinetic_energy_gradient()
  endif
  if (do_test_advect_density_field) then
    call run_advect_density_field()
  endif
  if (do_test_advect_theta_field) then
    call run_advect_theta_field()
  endif
  if (do_test_vorticity_advection) then
    call run_vorticity_advection()
  endif
  if (do_test_project_eos_pressure) then
    call run_project_eos_pressure()
  endif
  if (do_test_pressure_gradient_bd) then
    call run_pressure_gradient_bd()
  endif
  if (do_test_hydrostatic) then
    call run_hydrostatic()
  endif
  if (do_test_rk_alg) then
    call run_rk_alg()
  endif

  call finalise()

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program runge_kutta
