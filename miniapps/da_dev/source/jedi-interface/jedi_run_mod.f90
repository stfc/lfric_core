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

  use constants_mod,                 only : i_native

  implicit none

  private

type, public :: jedi_run_type
  private

contains

  !> Field initialiser.
  procedure, public :: initialise => jedi_run_initialiser

  !> Finalizer
  final             :: jedi_run_destructor

end type jedi_run_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> jedi_run constructor
!> @param [in] filename a character that contains the location of the namelist
!>             file
subroutine jedi_run_initialiser( self, filename )

  use mpi_mod,           only : initialise_comm
  use da_dev_driver_mod, only : initialise_lfric, initialise_lfric_comm

  implicit none

  class( jedi_run_type ), intent(inout) :: self
  character(len=*), intent(in)          :: filename

  integer(i_native) :: model_communicator
  integer(i_native) :: world_communicator

  ! JEDI will initialise MPI so calling it here to enforce that behaviour.
  ! It will be called outside the scope of the model interface.
  call initialise_comm(world_communicator)

  ! MPI has already been initialised so pass in the world communicator.
  call initialise_lfric_comm(model_communicator, world_communicator)

  ! initialise infrastructure
  call initialise_lfric(model_communicator, filename)

end subroutine jedi_run_initialiser

!> jedi_run finalizer
subroutine jedi_run_destructor(self)

  use da_dev_driver_mod, only : finalise_lfric

  implicit none

  type(jedi_run_type), intent(inout)    :: self

  call finalise_lfric()

end subroutine jedi_run_destructor

end module jedi_run_mod
