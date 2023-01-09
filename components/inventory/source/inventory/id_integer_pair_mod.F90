!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair integers with a unique identifier.
module id_integer_pair_mod

  use constants_mod,         only: i_def
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Integer Pair
  ! ========================================================================== !

  !> @brief An object pairing an integer with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_integer_pair_type

    private

    integer(kind=i_def) :: number_

  contains

    procedure, public :: initialise
    procedure, public :: get_integer

  end type id_integer_pair_type

contains

  !> @brief Initialises the id_integer_pair object
  !> @param[in] number   The number that will be stored in the paired object
  !> @param[in] id       The integer ID to pair with the integer
  subroutine initialise( self, number, id )

    implicit none

    class(id_integer_pair_type), intent(inout) :: self
    integer(kind=i_def),         intent(in)    :: number
    integer(kind=i_def),         intent(in)    :: id

    self%number_ = number
    call self%set_id(id)

  end subroutine initialise

  !> @brief Get the integer corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The integer
  function get_integer(self) result(number)

    implicit none

    class(id_integer_pair_type), target, intent(in) :: self
    integer(kind=i_def),                 pointer    :: number

    number => self%number_

  end function get_integer

end module id_integer_pair_mod
