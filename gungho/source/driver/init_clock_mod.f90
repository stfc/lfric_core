!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Prepare the clock component for use.
!>
module init_clock_mod

  use constants_mod,           only : r_second
  use calendar_mod,            only : calendar_type
  use clock_mod,               only : clock_type
  use log_mod,                 only : log_event,       &
                                      log_level_error, &
                                      log_scratch_space
  use step_calendar_mod,       only : step_calendar_type
  use time_config_mod,         only : calendar,          &
                                      calendar_timestep, &
                                      key_from_calendar, &
                                      timestep_end,      &
                                      timestep_start
  use timestepping_config_mod, only : dt, &
                                      spinup_period

  implicit none

  private
  public :: initialise_clock

contains

  !> Allocate a clock object and initialise from namelists.
  !>
  !> @param [out] clock Initialised clock. This should be NOT be allocated on
  !>                    entry.
  !>
  subroutine initialise_clock( clock )

    implicit none

    class(clock_type), intent(out), allocatable :: clock

    class(calendar_type), allocatable :: new_calendar
    integer :: status

    select case (calendar)
      case (calendar_timestep)
        allocate( step_calendar_type :: new_calendar, stat=status)
      case default
        write( log_scratch_space,                               &
               '("init_clock_mod: Unrecognised calendar: ",A)') &
            key_from_calendar( calendar )
        call log_event( log_scratch_space, log_level_error )
    end select ! case (calendar)
    if (status /= 0) then
      call log_event( 'init_clock_mod: Unable to allocate calendar', &
                      log_level_error )
    end if

    allocate( clock_type :: clock, stat=status)
    if (status /= 0) then
        call log_event( 'init_clock_mod: Unable to allocate clock', &
                        log_level_error )
    end if
    call clock%initialise( new_calendar,   &
                           timestep_start, &
                           timestep_end,   &
                           real(dt, kind=r_second), &
                           real(spinup_period, kind=r_second) )

  end subroutine initialise_clock

end module init_clock_mod
