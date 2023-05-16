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
  use driver_comm_mod,         only : init_comm, final_comm
  use driver_config_mod,       only : init_config, final_config
  use gravity_wave_mod,        only : program_name, &
                                      gravity_wave_required_namelists
  use gravity_wave_driver_mod, only : initialise, run, finalise
  use mpi_mod,                 only : global_mpi

  implicit none

  character(:), allocatable :: filename

  call init_comm( program_name )
  call get_initial_filename( filename )
  call init_config( filename, gravity_wave_required_namelists )
  deallocate( filename )
  call initialise( global_mpi )

  call run()

  call finalise()
  call final_config()
  call final_comm()

end program gravity_wave
