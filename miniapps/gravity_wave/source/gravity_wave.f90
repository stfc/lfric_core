!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page gravity_wave Gravity Wave miniapp
!> Test program for the automatic generation of boundary condition enforcement
!> by PSyclone.
!>
!> @brief Main program used to simulate the linear gravity waves equations.

program gravity_wave

  use cli_mod,                 only : get_initial_filename
  use gravity_wave_driver_mod, only : initialise, run, finalise
  use mod_wait,                only : init_wait
  use mpi_mod,                 only : initialise_comm, &
                                      finalise_comm
  use xios,                    only : xios_initialize, &
                                      xios_finalize

  implicit none

  character(*), parameter :: xios_id = "gravity_wave"

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999
  integer                   :: model_communicator = -999

  call get_initial_filename( filename )

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = model_communicator)

  call initialise( filename, model_communicator )
  deallocate( filename )

  call run()

  call finalise()

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program gravity_wave
