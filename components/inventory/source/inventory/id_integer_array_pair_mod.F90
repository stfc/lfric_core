!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair arrays of integers with a unique identifier.
module id_integer_array_pair_mod

  use constants_mod,         only: i_def
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Integer Array Pair
  ! ========================================================================== !

  !> @brief An object pairing an integer with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_integer_array_pair_type

    private

    integer(kind=i_def), allocatable :: numbers_(:)

  contains

    procedure, public :: initialise
    procedure, public :: get_integer_array

  end type id_integer_array_pair_type

contains

  !> @brief Initialises the id_integer_array_pair object
  !> @param[in] numbers  The array of numbers that will be stored in the paired object
  !> @param[in] id       The integer ID to pair with the integer
  subroutine initialise( self, numbers, id )

    implicit none

    class(id_integer_array_pair_type), intent(inout) :: self
    integer(kind=i_def),               intent(in)    :: numbers(:)
    integer(kind=i_def),               intent(in)    :: id

    integer(kind=i_def) :: i

    allocate(self%numbers_(size(numbers)))

    do i = 1, size(numbers)
      self%numbers_(i) = numbers(i)
    end do

    call self%set_id(id)

  end subroutine initialise

  !> @brief Get the integer corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The integer array
  function get_integer_array(self) result(numbers)

    implicit none

    class(id_integer_array_pair_type), target, intent(in) :: self
    integer(kind=i_def),                       pointer    :: numbers(:)

    numbers => self%numbers_

  end function get_integer_array

end module id_integer_array_pair_mod
