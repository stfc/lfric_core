!-------------------------------------------------------------------------------
!(c) Crown copyright 2020 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!>@brief Updates time-varying fields from linked list of time axis objects
module variable_fields_mod

  use constants_mod,                 only : r_def, r_second, str_def, i_def
  use clock_mod,                     only : clock_type
  use field_collection_mod,          only : field_collection_type
  use log_mod,                       only : log_event, &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO, &
                                            LOG_LEVEL_ERROR
  use lfric_xios_time_axis_mod,      only : time_axis_type
  use linked_list_mod,               only : linked_list_type, &
                                            linked_list_item_type
  implicit none

  private
  public :: init_variable_fields, &
            update_variable_fields

contains

  subroutine init_variable_fields(time_axis_list, clock, state)

    use xios, only: xios_date, xios_get_current_date, &
                    xios_date_get_day_of_year, xios_date_get_second_of_day

    implicit none

    type(linked_list_type),      intent(in)    :: time_axis_list
    class(clock_type),           intent(in)    :: clock
    type(field_collection_type), intent(inout) :: state

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()
    type(time_axis_type),        pointer :: time_axis => null()

    type(xios_date) :: datetime
    integer(i_def)  :: day_of_year
    real(r_def) :: time_sec, second_of_day

    ! Get the current time from xios. This will change when the clock is
    ! properly integrated.
    call xios_get_current_date(datetime)
    day_of_year   = int(xios_date_get_day_of_year(datetime), i_def)
    second_of_day = real(xios_date_get_second_of_day(datetime), r_def)
    time_sec = day_of_year * 24.0_r_def * 3600.0_r_def + second_of_day

    ! start at the head of the time_axis linked list
    loop => time_axis_list%get_head()

    do
      ! If list is empty or we're at the end of list, exit
      if ( .not. associated(loop) ) then
        exit
      end if

      ! Select the time_axis_type to get at the information in the list payload
      select type( list_item => loop%payload )
        type is (time_axis_type)
          time_axis => list_item

          ! Align time window and populate model data
          call time_axis%align(time_sec)
          call time_axis%update_fields()
          call time_axis%populate_model_fields(state, clock)

      end select

      loop => loop%next

    end do

    nullify(loop)
    nullify(time_axis)

  end subroutine init_variable_fields

  subroutine update_variable_fields(time_axis_list, clock, state)

    implicit none

    type(linked_list_type),      intent(in)    :: time_axis_list
    class(clock_type),           intent(in)    :: clock
    type(field_collection_type), intent(inout) :: state

    ! Pointer to linked list - used for looping through the list
    type(linked_list_item_type), pointer :: loop => null()
    type(time_axis_type),        pointer :: time_axis => null()

    real(r_def) :: time_window(2), time_sec

    ! Get time in seconds from clock
    ! Convert this to r_def because we will compare to time_window
    time_sec = real(clock%seconds_from_steps(clock%get_step()), r_def)

    ! Start at the head of the time_axis linked list
    loop => time_axis_list%get_head()

    do
      ! If list is empty or we're at the end of list, exit
      if ( .not. associated(loop) ) then
        exit
      end if

      ! Select the time_axis_type to get at the information in the list payload
      select type( list_item => loop%payload )
        type is (time_axis_type)
          time_axis => list_item

          write(log_scratch_space, '(A,A)') 'Updating variable fields from ', &
                                              trim(time_axis%get_name())
          call log_event(log_scratch_space, LOG_LEVEL_INFO)

          time_window = time_axis%get_time_window()

          ! If the time has moved to the next window, update the time axis
          if ( time_sec > time_window(2) ) then
            call time_axis%shift_forward()
            call time_axis%update_fields()
          end if

          ! If time is still out of bounds then something has gone wrong
          time_window = time_axis%get_time_window()
          if ( time_sec > time_window(2) ) then
            write(log_scratch_space, '(A,A)') 'Timestep larger than time axis resolution: ', &
                                              trim(time_axis%get_name())
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          end if

          ! Populate the model fields from the time axis data
          if ( mod(clock%get_step(), time_axis%get_update_frequency()) == 0 ) then
            call time_axis%populate_model_fields(state, clock)
          end if

      end select

      loop => loop%next

    end do

    nullify(loop)
    nullify(time_axis)

  end subroutine update_variable_fields

end module variable_fields_mod
