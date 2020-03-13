!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page transport Transport miniapp
!> Program file for running transport miniapp. Subroutine calls include initialise_transport(),
!> run_transport() and finalise_transport().
program transport

  use cli_mod,              only: get_initial_filename
  use mod_wait,             only: init_wait
  use mpi_mod,              only: finalise_comm, &
                                  initialise_comm
  use transport_driver_mod, only: initialise_transport, &
                                  run_transport,        &
                                  finalise_transport
  use xios,                 only: xios_finalize, &
                                  xios_initialize

  implicit none

  character(len=*), parameter :: xios_id = "transport"

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999
  integer                   :: model_communicator = -999

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  ! Initialise XIOS and get back the split mpi communicator
  call init_wait()
  call xios_initialize(xios_id, return_comm = model_communicator)

  call get_initial_filename( filename )
  call initialise_transport( filename, model_communicator )
  deallocate( filename )

  call run_transport()

  call finalise_transport()

  ! Finalise XIOS
  call xios_finalize()

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program transport
