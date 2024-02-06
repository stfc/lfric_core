!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Container for time axes.
!>
module gungho_time_axes_mod

  use key_value_collection_mod,     only : key_value_collection_type
  use key_value_mod,                only : abstract_value_type
  use linked_list_mod,              only : linked_list_type

  implicit none

  private
  public :: get_time_axes_from_collection

  !> Collection of time axes.
  !>
  type, extends(abstract_value_type), public :: gungho_time_axes_type

    private

    !> Time varying ancillaries time axis.
    type(linked_list_type), public :: ancil_times_list

    !> Time varying LBC time axis.
    type(linked_list_type), public :: lbc_times_list

    !> Time varying linearisation state time axis.
    !>
    !> @todo Is this part of the linear model?
    !>
    type(linked_list_type), public :: ls_times_list

  end type gungho_time_axes_type

contains

  !-----------------------------------------------------------------------------
  ! Non-type-bound helper function
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> @brief Helper function to extract a concrete time axes object from a
  !>        key-value collection
  !> @param[in] collection The key-value collection to extract from
  !> @param[in] name       The name of the time axes object to extract
  !> @return    fs_chain   The requested time_axes object
  function get_time_axes_from_collection(collection, name) result(time_axes)

  implicit none

    type(key_value_collection_type), intent(in) :: collection
    character(*),                    intent(in) :: name

    type(gungho_time_axes_type), pointer    :: time_axes

    class(abstract_value_type), pointer :: abstract_value

    call collection%get_value(trim(name), abstract_value)
    select type(abstract_value)
      type is (gungho_time_axes_type)
      time_axes => abstract_value
    end select

  end function get_time_axes_from_collection

end module gungho_time_axes_mod
