!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page gung_ho GungHo Program
!> This is a code that uses the LFRic infrastructure to build a model that
!> just includes the GungHo dynamical core.

!> @brief Main program used to illustrate gungho functionality.

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program gungho

  use cli_mod,           only : get_initial_filename
  use gungho_driver_mod, only : initialise, run, finalise
  use mod_wait,          only : init_wait
  use mpi_mod,           only : initialise_comm, &
                                finalise_comm
  use xios,              only : xios_initialize, &
                                xios_finalize

  implicit none

  character(*), parameter :: xios_id = "gungho"

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

end program gungho
