!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> A module containing utility functions to contain events to be attached to
!> the ticking of a model clock
!>
module event_actor_mod

  use constants_mod,        only : str_def
  use linked_list_data_mod, only : linked_list_data_type
  use log_mod,              only : log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private

  ! Abstract class used for timestep event interface
  type, extends(linked_list_data_type), public, abstract :: event_actor_type
    private
    character(str_def) :: event_name
    logical :: active = .false.
    logical :: constructed = .false.
  contains
    procedure, public :: init_event_actor
    procedure, public :: get_event_name
    procedure, public :: is_active
    procedure, public :: set_active
  end type event_actor_type

contains

  subroutine init_event_actor(this, name)
    implicit none
    class(event_actor_type), intent(inout) :: this
    character(*), intent(in) :: name

    if(this%constructed) then
      write(log_scratch_space, '(A)') trim(name) // &
                                      " event actor type already initialised"
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    this%event_name = name // "_event"
    this%active = .true.
    this%constructed = .true.
  end subroutine init_event_actor

  !> Returns the name of the event
  function get_event_name(this) result(name)

    implicit none

    class(event_actor_type), intent(in) :: this
    character(str_def) :: name

    name = this%event_name

  end function get_event_name

  !> Returns true if the event actor can be used
  function is_active(this) result(l_active)

    implicit none

    class(event_actor_type), intent(inout) :: this
    logical                                :: l_active

    l_active = this%active

  end function is_active

  !> Sets active state of event actor
  !> @param [in] active Set active state of event actor
  subroutine set_active(this, active)

    implicit none

    class(event_actor_type), intent(inout) :: this
    logical, intent(in) :: active

    this%active = active

  end subroutine

end module event_actor_mod
