!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page linear Linear model Program
!> This is a code that uses the LFRic infrastructure to build a model that
!> is the tangent linear/ perturbation forecast for gungho.

!> @brief Main program used to illustrate linear model functionality.

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program linear_model

  use cli_mod,               only : get_initial_filename
  use driver_comm_mod,       only : init_comm, final_comm
  use driver_config_mod,     only : init_config, final_config
  use driver_log_mod,        only : init_logger, final_logger
  use gungho_mod,            only : gungho_required_namelists
  use gungho_modeldb_mod,    only : modeldb_type
  use linear_driver_mod,     only : initialise, run, finalise
  use mpi_mod,               only : global_mpi

  implicit none

  ! Model run working data set
  type (modeldb_type) :: modeldb

  character(*), parameter :: application_name = "linear_model"

  character(:), allocatable :: filename

  modeldb%mpi => global_mpi

  call init_comm( application_name )
  call get_initial_filename( filename )
  call init_config( filename, gungho_required_namelists )
  deallocate( filename )
  call init_logger( modeldb%mpi%get_comm(), application_name)

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%model_data%depository%initialise(name='depository', table_len=100)
  call modeldb%model_data%prognostic_fields%initialise(name="prognostics", table_len=100)
  call modeldb%model_data%diagnostic_fields%initialise(name="diagnostics", table_len=100)

  call initialise( application_name, modeldb )
  call run( application_name, modeldb )
  call finalise( application_name, modeldb )

  call final_logger( application_name )
  call final_config()
  call final_comm()

end program linear_model
