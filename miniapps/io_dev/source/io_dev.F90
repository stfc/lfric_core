!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page io_dev io_dev Miniapp
!> Test program for the XIOS IO implementation.

!> @brief Main program used to test XIOS setup and output of a field.

program io_dev

  use cli_mod,           only: get_initial_filename
  use io_dev_driver_mod, only: initialise, run, finalise
  use mod_wait,          only: init_wait
  use mpi_mod,           only: finalise_comm, &
                               initialise_comm
  use xios,              only: xios_finalize, &
                               xios_initialize

  implicit none

  character(*), parameter :: xios_id = 'io_dev'

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999
  integer                   :: model_communicator = -999

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize( xios_id, return_comm=model_communicator )

  call get_initial_filename( filename )
  call initialise ( filename, model_communicator )
  deallocate ( filename )

  call run()

  call finalise()

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program io_dev
