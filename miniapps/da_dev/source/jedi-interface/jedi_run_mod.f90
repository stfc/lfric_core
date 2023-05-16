!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing an class that handles LFRic initialisation.
!>
!> @details This class handles the initialisation and finalisation of LFRic
!
module jedi_run_mod

  use constants_mod, only : i_native, str_def
  use mpi_mod,       only : mpi_type

  implicit none

  private

type, public :: jedi_run_type
  private
  character(str_def) :: jedi_run_name

contains

  !> Field initialiser.
  procedure, public :: initialise

  !> Finalizer
  final             :: jedi_run_destructor

end type jedi_run_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_run_type
!>
!> @param [in] filename a character that contains the location of the namelist
!>             file
subroutine initialise( self, program_name, mpi )

  use mpi_mod,           only : create_comm
  use da_dev_driver_mod, only : initialise_lfric

  implicit none

  class( jedi_run_type ), intent(inout) :: self
  character(len=*), intent(in)          :: program_name
  class(mpi_type),  intent(inout)       :: mpi

  self%jedi_run_name = program_name

  ! initialise infrastructure
  call initialise_lfric( program_name, mpi )

end subroutine initialise

!> @brief    Finalizer for jedi_run_type
!>
subroutine jedi_run_destructor(self)

  use da_dev_driver_mod, only : finalise_lfric

  implicit none

  type(jedi_run_type), intent(inout)    :: self

  call finalise_lfric( trim(self%jedi_run_name) )

end subroutine jedi_run_destructor

end module jedi_run_mod
