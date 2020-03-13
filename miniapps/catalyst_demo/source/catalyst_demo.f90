!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page catalyst ParaView Catalyst visualisation mini app
!> Test program for demonstrating in-situ visualisation with ParaView
!> Catalyst.
!>
!> @brief Main program used to run an example for in-situ visualisation.

program catalyst_demo

  use catalyst_demo_driver_mod, only: initialise, run, finalise
  use cli_mod,                  only: get_initial_filename
  use mod_wait,                 only: init_wait
  use mpi_mod,                  only: finalise_comm, &
                                      initialise_comm

  implicit none

  character(*), parameter :: xios_id = 'catalyst_demo'

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999
  integer                   :: model_communicator = -999

  call get_initial_filename( filename )

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = model_communicator)

  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program catalyst_demo
