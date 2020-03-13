!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp cell_locator program

!> @brief Finds the cell index that contains a target point

!> @details Calls init, run and finalise routines from a driver module

program cell_locator

  use cell_locator_driver_mod, only : initialise, run, finalise
  use cli_mod,                 only : get_initial_filename
  use mpi_mod,                 only : finalise_comm, &
                                      initialise_comm

  implicit none

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999

  ! Initialise mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  call get_initial_filename( filename )
  call initialise( filename, world_communicator )
  deallocate( filename )

  call run()

  call finalise()

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program cell_locator
