!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair field arrays with a unique identifier.
module id_r32_field_array_pair_mod

  use constants_mod,         only: i_def
  use field_r32_mod,         only: field_r32_type
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Field Array Pair
  ! ========================================================================== !

  !> @brief An object pairing a field array with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_r32_field_array_pair_type

    private

    type(field_r32_type), allocatable :: field_array_(:)

  contains

    procedure, public :: initialise
    procedure, public :: get_field_array

  end type id_r32_field_array_pair_type

contains

  !> @brief Initialises the id_r32_field_array_pair object
  !> @param[in] field_array   The fields to be stored in the paired object
  !> @param[in] id            The integer ID to pair with the field_array
  subroutine initialise( self, field_array, id )

    implicit none

    class(id_r32_field_array_pair_type), intent(inout) :: self
    type(field_r32_type),                intent(in)    :: field_array(:)
    integer(kind=i_def),                 intent(in)    :: id

    integer(kind=i_def) :: i

    allocate(self%field_array_(size(field_array)))

    do i = 1, size(field_array)
      call self%field_array_(i)%initialise( vector_space=field_array(i)%get_function_space() )
      call field_array(i)%copy_field(self%field_array_(i))
    end do

    call self%set_id(id)

  end subroutine initialise

  !> @brief Get the field_array corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The field_array
  function get_field_array(self) result(field_array)

    implicit none

    class(id_r32_field_array_pair_type), target, intent(in) :: self
    type(field_r32_type),                        pointer    :: field_array(:)

    field_array => self%field_array_

  end function get_field_array

end module id_r32_field_array_pair_mod
