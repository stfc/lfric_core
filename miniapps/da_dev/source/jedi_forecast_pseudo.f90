!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page jedi_forecast_pseudo program

!> @brief Main program for running pseudo model forecast with jedi emulator
!>        objects.

!> @details Setup and run a pseudo model forecast using the jedi emulator
!>          objects. The jedi objects are constructed via an initialiser call
!>          and the forecast is handled by the model object.
!>
program jedi_forecast_pseudo

  use constants_mod,     only : i_def
  use da_dev_mod,        only : da_dev_required_namelists
  use da_dev_driver_mod, only : finalise_model
  use driver_comm_mod,   only : init_comm, final_comm
  use driver_config_mod, only : init_config, final_config
  use mpi_mod,           only : global_mpi

  ! Data types and methods to get/store configurations
  use jedi_state_config_mod,        only : jedi_state_config_type
  use jedi_pseudo_model_config_mod, only : jedi_pseudo_model_config_type
  use cli_mod,                      only : get_initial_filename

  ! Jedi emulator objects
  use jedi_run_mod,          only : jedi_run_type
  use jedi_geometry_mod,     only : jedi_geometry_type
  use jedi_state_mod,        only : jedi_state_type
  use jedi_pseudo_model_mod, only : jedi_pseudo_model_type


  implicit none

  ! Jedi objects
  type(jedi_geometry_type)     :: jedi_geometry
  type(jedi_state_type)        :: jedi_state
  type(jedi_pseudo_model_type) :: jedi_model
  type(jedi_run_type)          :: jedi_run

  ! Emulator configs
  type(jedi_state_config_type)        :: jedi_state_config
  type(jedi_pseudo_model_config_type) :: jedi_pseudo_model_config
  integer( kind=i_def )               :: date_time_duration
  character(:), allocatable           :: filename

  character(*), parameter      :: program_name = "jedi_forecast_pseudo"

  ! Infrastructure config
  call init_comm( program_name )
  call get_initial_filename( filename )
  call init_config( filename, da_dev_required_namelists )
  deallocate( filename )

  ! Run object
  ! Handles initialization and finalization of required infrastructure
  call jedi_run%initialise( program_name, global_mpi )

  ! Config for the jedi emulator objects
  ! State config
  call jedi_state_config%initialise( use_nl_model = .false. )

  ! Model config
  call jedi_pseudo_model_config%initialise()

  ! Forecast config - duration of forecast / seconds
  date_time_duration = 5_i_def

  ! Geometry
  call jedi_geometry%initialise()

  ! State
  call jedi_state%initialise( jedi_geometry, jedi_state_config )

  ! Model
  call jedi_model%initialise( jedi_pseudo_model_config )

  ! Run app via model class
  call jedi_model%forecast( jedi_state, date_time_duration )

  ! To provide KGO
  call finalise_model( program_name, jedi_state%io_collection )

  call final_config()
  call final_comm()

end program jedi_forecast_pseudo
