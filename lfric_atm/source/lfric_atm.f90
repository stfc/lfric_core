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
  use gungho_driver_mod, only : initialise, run, finalise

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

end program lfric_atm
