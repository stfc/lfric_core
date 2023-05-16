!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp diagnostics program

!> @brief Program used as proof of concept for testing core infrastructure. Simple miniapp for easy cannibalising.

!> @details Calls init, run and finalise routines from a driver module

program diagnostics

    use cli_mod,                       only : get_initial_filename
    use diagnostics_configuration_mod, only : program_name
    use driver_comm_mod,               only : init_comm, final_comm
    use driver_config_mod,             only : init_config, final_config
    use diagnostics_configuration_mod, only : required_namelists
    use diagnostics_driver_mod,        only : initialise, run, finalise
    use mpi_mod,                       only : global_mpi

    implicit none

    character(:), allocatable :: filename

    call init_comm( program_name )
    call get_initial_filename( filename )
    call init_config( filename, required_namelists )
    deallocate( filename )
    call initialise( global_mpi )

    call run()

    call finalise()
    call final_config()
    call final_comm()

end program diagnostics
