!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class than handles the pseudo model
!>
!> @details This module includes a class that handles the pseudo model where
!>          the time-stepping uses the read_file state method. An included
!>          forecast application uses the model init, step and final to run a
!>          forecast.
!>
module jedi_pseudo_model_mod

  use jedi_state_mod,                only : jedi_state_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use constants_mod,                 only : i_def, str_def

  implicit none

  private

type, public :: jedi_pseudo_model_type
  private

  !> A date_time list to read
  integer( kind=i_def ), allocatable :: date_time_states(:)
  !> index for the current state in date_time_states
  integer( kind=i_def )               :: current_state
  !> the number of states
  integer( kind=i_def )               :: n_states
  !> the file prefix for reading
  character(len=str_def)              :: file_prefix

contains

  !> Model initialiser.
  procedure, public :: initialise => jedi_pseudo_model_initialiser

  !> Methods
  procedure, private :: model_init
  procedure, private :: model_step
  procedure, private :: model_final

  !> Public method
  procedure, public :: forecast

  !> Finalizer
  final             :: jedi_pseudo_model_destructor

end type jedi_pseudo_model_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> jedi_pseudo_model constructor
!> @param [in] config configuration used to setup the model class
subroutine jedi_pseudo_model_initialiser( self, config)

  use jedi_pseudo_model_config_mod, only : jedi_pseudo_model_config_type

  implicit none

  class( jedi_pseudo_model_type ), intent(inout)  :: self
  type(jedi_pseudo_model_config_type), intent(in) :: config

  ! Setup the pseudo model
  self%current_state = 1
  self%n_states = size(config%date_time_states, dim=1)
  allocate(self%date_time_states(self%n_states))
  self%date_time_states = config%date_time_states
  self%file_prefix = config%read_file_prefix

end subroutine jedi_pseudo_model_initialiser

!> Model initialise
!> @param [inout] state state object to be used in the model initialise
subroutine model_init(self, state)

  implicit none

  class( jedi_pseudo_model_type ), intent(in) :: self
  type( jedi_state_type ),  intent(inout)     :: state

end subroutine model_init

!> Model step
!> @param [inout] state state object to propagated
subroutine model_step(self, state)

  implicit none
  class( jedi_pseudo_model_type ), intent(inout) :: self
  type( jedi_state_type ),  intent(inout)        :: state

  ! Local
  integer :: date_time

  if (self%current_state>self%n_states) then
    write(log_scratch_space, '(A)') "self%current_state>self%n_states."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  date_time = self%date_time_states(self%current_state)
  call state%read_file(date_time, self%file_prefix)

  ! Iterate the current state
  self%current_state = self%current_state + 1

end subroutine model_step

!> Model finalise
!> @param [inout] state state object to be used in the model finalialise
subroutine model_final(self, state)

  implicit none

  class( jedi_pseudo_model_type ), intent(in) :: self
  type( jedi_state_type ),  intent(inout) :: state

end subroutine model_final

!> Finalizer
subroutine jedi_pseudo_model_destructor(self)

  implicit none

  type(jedi_pseudo_model_type), intent(inout) :: self

  if (allocated(self%date_time_states)) deallocate(self%date_time_states)

end subroutine jedi_pseudo_model_destructor

!-------------------------------------------------------------------------------
! OOPS defined forecast method
!-------------------------------------------------------------------------------

!> Model forecast
!> @param [inout] state the state object to propagate
!> @param [in] date_time_duration the duration of the forecast
subroutine forecast(self, state, date_time_duration)

  implicit none

  class( jedi_pseudo_model_type ), intent(inout) :: self
  type( jedi_state_type ),  intent(inout)        :: state
  integer, intent(in)                            :: date_time_duration

  ! Local
  integer :: date_time_end

  ! end time
  date_time_end = state%valid_time() + date_time_duration

  call self%model_init(state)

  ! initialize the post processor
  ! post.initialize(state, ...);
  ! In a H(X) the following call would be used to do spatial interpolations
  ! on the Atlas/dummy field data
  ! post%process(state)

  ! loop until date_time_end
  do while (state%valid_time() < date_time_end)
    call self%model_step(state)
    ! post%process(state)
  enddo

  ! post.finalize(state);
  call self%model_final(state)

end subroutine forecast

end module jedi_pseudo_model_mod
