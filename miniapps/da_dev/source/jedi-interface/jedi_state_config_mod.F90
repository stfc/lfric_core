!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a faux config for the mock jedi state.
!>
!> @details A class is defined to hold the configuration data required to
!>          construct a mock jedi state. An initialiser is included and
!>          this is currently hard coded. In JEDI this information would be
!>          stored in a yaml file and eckit is used to parse and store a
!>          configuration object.
!>
module jedi_state_config_mod

  use constants_mod,      only : i_def, str_def, l_def
  use fs_continuity_mod,  only : W3, Wtheta

  implicit none

  private

type, public :: jedi_state_config_type

  !> the name of the state variables
  character( len=str_def ), allocatable :: state_variables(:)
  !> function space for the variables
  integer( kind=i_def ), allocatable :: variable_function_space(:)
   !> logical that defines if the field is 2D
  logical( kind=l_def ), allocatable :: variable_is_2d(:)

  ! here we have date_time as an integer. It will actually be an object or string that stores time
  ! to be read. initilly stored in the configaurtion file:
  !
  ! date: '2018-04-14T21:00:00Z'
  ! mpas define the formating for this as:
  ! dateTimeString = '$Y-$M-$D_$h:$m:$s'
  integer( kind=i_def )              :: date_time
  !> file prefix for read
  character(len=str_def)             :: read_file_prefix
  !> Defines if the model_data is ...
  logical( kind=l_def )              :: use_full_model

contains

  !> jedi_state_config initialiser.
  procedure :: initialise

  !> jedi_state_config finalizer
  final     :: jedi_state_config_destructor

end type jedi_state_config_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> jedi_state_config initialiser
!> @param [in] use_full_model set true if full model is used
subroutine initialise( self, use_full_model )

  implicit none

  class( jedi_state_config_type ), intent(inout) :: self
  logical( kind=l_def ), intent(in)              :: use_full_model
  ! Local
  integer :: nvars

  ! configuration inputs
  self%read_file_prefix="read_"
  nvars=2
  allocate(self%state_variables(nvars), self%variable_function_space(nvars), self%variable_is_2d(nvars))
  !
  self%state_variables(1) = "theta"
  self%state_variables(2) = "rho"
  !self%state_variables(3) = "u"
  !
  self%variable_function_space(1) = Wtheta
  self%variable_function_space(2) = W3
  !self%variable_function_space(3) = W2
  !
  self%variable_is_2d(1) = .false.
  self%variable_is_2d(2) = .false.
  !self%variable_is_2d(3) = .false.
  self%date_time = 0

  self%use_full_model = use_full_model

end subroutine initialise

!> jedi_state_config finalizer
subroutine jedi_state_config_destructor(self)!

  implicit none

  type(jedi_state_config_type), intent(inout)    :: self

  if ( allocated(self%state_variables) ) deallocate(self%state_variables)
  if ( allocated(self%variable_function_space) ) deallocate(self%variable_function_space)
  if ( allocated(self%variable_is_2d) ) deallocate(self%variable_is_2d)

end subroutine jedi_state_config_destructor

end module jedi_state_config_mod
