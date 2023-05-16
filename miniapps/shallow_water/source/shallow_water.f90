!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page shallow_water Shallow water equations miniapp
!> This is code that uses the LFRic infrastructure to build a shallow water
!> model that includes some of the GungHo routines.
!>
!> @brief Main program used to simulate shallow water equations.
!>
!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program shallow_water

  use cli_mod,                  only: get_initial_filename
  use driver_comm_mod,          only: init_comm, final_comm
  use driver_config_mod,        only: init_config, final_config
  use mpi_mod,                  only: global_mpi
  use shallow_water_mod,        only: program_name, &
                                      shallow_water_required_namelists
  use shallow_water_driver_mod, only: initialise, &
                                      run,        &
                                      finalise

  implicit none

  character(:), allocatable :: filename

  call init_comm( program_name )
  call get_initial_filename( filename )
  call init_config( filename, shallow_water_required_namelists )
  deallocate( filename )

  call initialise( global_mpi )

  call run()

  call finalise()
  call final_config()
  call final_comm()

end program shallow_water
