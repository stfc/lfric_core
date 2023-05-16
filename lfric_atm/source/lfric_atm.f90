!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> This is a code that uses the LFRic infrastructure to build a model that
!> includes the GungHo dynamical core and physics parametrisation schemes
!> that are currently provided through the use of unified model code.

!> @brief Main program used to illustrate an atmospheric model built using
!>        LFRic infrastructure

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the atmospheric model.

program lfric_atm

  use cli_mod,           only : get_initial_filename
  use driver_comm_mod,   only : init_comm, final_comm
  use driver_config_mod, only : init_config, final_config
  use gungho_mod,        only : gungho_required_namelists
  use gungho_driver_mod, only : initialise, run, finalise
  use mpi_mod,           only : global_mpi

  implicit none

  character(*), parameter :: application_name = "lfric_atm"

  character(:), allocatable :: filename

  call init_comm( application_name )
  call get_initial_filename( filename )
  call init_config( filename, gungho_required_namelists )
  deallocate( filename )
  call initialise( application_name, global_mpi )

  call run( application_name )

  call finalise( application_name )
  call final_config()
  call final_comm()

end program lfric_atm
