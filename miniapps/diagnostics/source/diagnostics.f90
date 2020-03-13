!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp diagnostics program

!> @brief Program used as proof of concept for testing core infrastructure. Simple miniapp for easy cannibalising.

!> @details Calls init, run and finalise routines from a driver module

program diagnostics

    use cli_mod,                only : get_initial_filename
    use diagnostics_driver_mod, only : initialise, run, finalise
    use mod_wait,               only : init_wait
    use mpi_mod,                only : finalise_comm, &
                                       initialise_comm
    use xios,                   only : xios_finalize, &
                                       xios_initialize

    use log_mod, only : log_event, LOG_LEVEL_INFO

    implicit none

    character(*), parameter :: xios_id = "diagnostics"

    character(:), allocatable :: filename
    integer :: model_communicator = -999
    integer :: world_communicator = -999

    call get_initial_filename(filename)

    ! Initialse mpi and create the default communicator: mpi_comm_world
    call initialise_comm( world_communicator )

    ! Initialise XIOS and get back the split mpi communicator
    call init_wait()
    call xios_initialize( xios_id, return_comm=model_communicator )

    call initialise( filename, model_communicator )
    deallocate(filename)

    !do it!
    call run()

    call finalise()

    ! Finalise XIOS
    call xios_finalize()

    ! Finalise mpi and release the communicator
    call finalise_comm()

end program diagnostics
