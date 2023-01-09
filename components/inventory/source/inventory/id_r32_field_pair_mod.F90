!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair r32 fields with a unique identifier.
module id_r32_field_pair_mod

  use constants_mod,         only: i_def
  use field_r32_mod,         only: field_r32_type
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Field Pair
  ! ========================================================================== !

  !> @brief An object pairing an operator_type with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_r32_field_pair_type

    private

    type(field_r32_type) :: field_

  contains

    procedure, public :: initialise
    procedure, public :: get_field

  end type id_r32_field_pair_type

contains

  !> @brief Initialises the id_r32_field_pair object
  !> @param[in] field    The field that will be stored in the paired object
  !> @param[in] id       The integer ID to pair with the field
  subroutine initialise( self, field, id )

    implicit none

    class(id_r32_field_pair_type), intent(inout) :: self
    type(field_r32_type),          intent(in)    :: field
    integer(kind=i_def),           intent(in)    :: id

    call self%field_%initialise( vector_space=field%get_function_space() )
    call field%copy_field(self%field_)
    call self%set_id(id)

  end subroutine initialise

  !> @brief Get the field corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The field
  function get_field(self) result(field)

    implicit none

    class(id_r32_field_pair_type), target, intent(in) :: self
    type(field_r32_type),                  pointer    :: field

    field => self%field_

  end function get_field

end module id_r32_field_pair_mod
