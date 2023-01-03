!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a class than handles the non-linear model
!>
!> @details This module includes a class that handles the non-linear models
!>          time stepping. An included forecast application uses the model
!>          init, step and final to run a forecast.
!>
module jedi_model_mod

  use jedi_state_mod,                only : jedi_state_type
  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use constants_mod,                 only : i_def

  implicit none

  private

type, public :: jedi_model_type
  private

  !> The data map between external field data and LFRic fields
  integer( kind = i_def )            :: date_time_duration_dt

contains

  !> Model initialiser.
  procedure, public :: initialise => jedi_model_initialiser

  !> methods
  procedure, private :: model_init
  procedure, private :: model_step
  procedure, private :: model_final

  !> public method
  procedure, public :: forecast

  !> Finalizer
  final             :: jedi_model_destructor

end type jedi_model_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> jedi_model constructor
!> @param [in] date_time_duration_dt time step duration
subroutine jedi_model_initialiser( self, date_time_duration_dt )

  implicit none

  class( jedi_model_type ), intent(inout)   :: self
  integer( kind=i_def ), intent(in)         :: date_time_duration_dt

  self%date_time_duration_dt = date_time_duration_dt

end subroutine jedi_model_initialiser

!> Model initialise
!> @param [in] state state object to be used in the model initialise
subroutine model_init(self, state)

  implicit none

  class( jedi_model_type ), target, intent(in) :: self
  type( jedi_state_type ),  intent(inout)      :: state

end subroutine model_init

!> Model step
!> @param [in] state state object to propagated by one time-step
subroutine model_step(self, state)

  use da_dev_driver_mod,  only: model_clock
  use da_dev_driver_mod,  only: step_lfric

  implicit none

  class( jedi_model_type ), target, intent(in) :: self
  type( jedi_state_type ),  intent(inout)      :: state

  ! Local
  logical :: clock_stopped

  ! check the clock
  clock_stopped=.not.model_clock%tick()
  !! if the clock has finished then it will just get the
  !! data at the end of the file - this prevents that
  if (clock_stopped) then
    write(log_scratch_space, '(A)') "Model::model_step::The LFRic clock has stopped."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  call step_lfric(state%model_data)

  ! Link the state fields with the model data
  call state%setup_external()
  call state%from_model_data()

  ! update the state time
  call state%update_time(self%date_time_duration_dt)

end subroutine model_step

!> Model finalise
!> @param [in] state state object to be used in the model finalialise
subroutine model_final(self, state)

  implicit none

  class( jedi_model_type ), target, intent(in) :: self
  type( jedi_state_type ),  intent(inout)      :: state

end subroutine model_final

!> Finalizer
subroutine jedi_model_destructor(self)

  implicit none

  type(jedi_model_type), intent(inout) :: self

end subroutine jedi_model_destructor

!-------------------------------------------------------------------------------
! OOPS defined forecast method
!-------------------------------------------------------------------------------

!> Model forecast
!> @param [inout] state the state object to propagate
!> @param [in] date_time_duration the duration of the forecast
subroutine forecast(self, state, date_time_duration)

  implicit none

  class( jedi_model_type ), target, intent(in) :: self
  type( jedi_state_type ),  intent(inout)      :: state
  integer, intent(in)                          :: date_time_duration

  ! Local
  integer :: date_time_end

  ! End time
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

end module jedi_model_mod
