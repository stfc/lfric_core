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

program lfric_coupled

  use cli_mod,            only : get_initial_filename
  use gungho_driver_mod,  only : initialise, run, finalise
  use mod_wait,           only : init_wait
  use mpi_mod,            only : initialise_comm, &
                                 finalise_comm
  use xios,               only : xios_initialize, &
                                 xios_finalize
  use derived_config_mod, only : l_esm_couple
  use coupler_mod,        only : cpl_initialize, cpl_finalize

  implicit none

  character(*), parameter :: xios_id   = "lfric_atmosphere"

  character(:), allocatable :: filename
  integer                   :: world_communicator = -999
  integer                   :: model_communicator = -999

  call get_initial_filename( filename )

  ! Initialse mpi and create the default communicator: mpi_comm_world
  call initialise_comm( world_communicator )

  ! Initialise XIOS and get back the split mpi communicator
  if(l_esm_couple) then
   !mpi already initialized
       call cpl_initialize( model_communicator )
       call xios_initialize(xios_id, local_comm = model_communicator)
  else
       call init_wait()
       call xios_initialize(xios_id, return_comm = model_communicator)
  endif


  call initialise( filename, model_communicator )
  deallocate( filename )

  call run()

  call finalise()

  ! Finalise XIOS
  call xios_finalize()

  if(l_esm_couple) then
     call cpl_finalize()
  endif

  ! Finalise mpi and release the communicator
  call finalise_comm()

end program lfric_coupled
