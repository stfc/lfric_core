!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a faux config for the mock pseudo model.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct the mock pseudo model. An initialiser is included and
!>          this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_pseudo_model_config_mod

  use constants_mod,      only : i_def, str_def

  implicit none

  private

type, public :: jedi_pseudo_model_config_type

  ! here we have date_time as an integer. It will actually be an object or string that stores time
  ! to be read. initilly stored in the configaurtion file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'

  !> List of the dates to be read
  integer( kind=i_def ), allocatable :: date_time_states(:)

  !> file prefix for read
  character(len=str_def)             :: read_file_prefix

contains

  !> jedi_pseudo_model initialiser.
  procedure :: initialise

  !> Finalizer
  final     :: jedi_pseudo_model_config_destructor

end type jedi_pseudo_model_config_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> jedi_pseudo_model_config initialiser
subroutine initialise( self )

  implicit none

  class( jedi_pseudo_model_config_type ), intent(inout) :: self
  ! Local
  integer( kind=i_def ) :: date_time_entries

  date_time_entries = 9
  allocate(self%date_time_states(date_time_entries))
  self%date_time_states = (/1,2,3,4,5,6,7,8,9/)
  !self%date_time_states = (/1,5,9/)

  self%read_file_prefix="read_"

end subroutine initialise

!> jedi_pseudo_model_config finalizer
subroutine jedi_pseudo_model_config_destructor(self)!

  implicit none

  type(jedi_pseudo_model_config_type), intent(inout)    :: self

  if ( allocated(self%date_time_states) ) deallocate(self%date_time_states)

end subroutine jedi_pseudo_model_config_destructor

end module jedi_pseudo_model_config_mod
