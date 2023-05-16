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

program gungho_model

  use cli_mod,           only : get_initial_filename
  use driver_comm_mod,   only : init_comm, final_comm
  use driver_config_mod, only : init_config, final_config
  use gungho_mod,        only : gungho_required_namelists
  use gungho_driver_mod, only : initialise, run, finalise
  use mpi_mod,           only : global_mpi

  implicit none

  character(*), parameter :: application_name = "gungho_model"

  character(:), allocatable :: filename

  call init_comm( application_name )
  call get_initial_filename( filename )
  call init_config( filename, gungho_required_namelists )

  call initialise( application_name, global_mpi )

  call run( application_name )

  call finalise( application_name )
  call final_config()
  call final_comm()

end program gungho_model
