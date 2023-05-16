!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp da_dev program

!> @brief Main program for running da_dev independently.

!> @details Calls init, run and finalise routines from a driver module

program da_dev

  use cli_mod,                   only : get_initial_filename
  use da_dev_mod,                only : da_dev_required_namelists
  use da_dev_driver_mod,         only : initialise_lfric, run, finalise_lfric, &
                                        initialise_model, finalise_model
  use driver_comm_mod,           only : init_comm, final_comm
  use driver_config_mod,         only : init_config, final_config
  use driver_model_data_mod,     only : model_data_type
  use constants_mod,             only : i_native
  use mpi_mod,                   only : global_mpi

  implicit none

  type(model_data_type) :: model_data

  character(*), parameter :: program_name = "da_dev"

  character(:), allocatable :: filename

  call init_comm( program_name )
  call get_initial_filename( filename )
  call init_config( filename, da_dev_required_namelists )
  deallocate( filename )
  call initialise_lfric( program_name, global_mpi )

  call initialise_model(global_mpi, model_data)

  call run(program_name, model_data)

  call finalise_model(program_name, model_data%depository)

  call finalise_lfric(program_name)
  call final_config()
  call final_comm()

end program da_dev
