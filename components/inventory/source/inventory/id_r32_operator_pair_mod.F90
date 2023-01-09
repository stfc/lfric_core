!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Defines an object to pair operators with a unique identifier.
module id_r32_operator_pair_mod

  use constants_mod,         only: i_def
  use operator_r32_mod,      only: operator_r32_type
  use id_abstract_pair_mod,  only: id_abstract_pair_type

  implicit none

  private

  ! ========================================================================== !
  ! ID-Operator Pair
  ! ========================================================================== !

  !> @brief An object pairing a field with a unique identifier
  !>
  type, public, extends(id_abstract_pair_type) :: id_r32_operator_pair_type

    private

    type(operator_r32_type) :: operator_

  contains

    procedure, public :: initialise
    procedure, public :: get_operator

  end type id_r32_operator_pair_type

contains

  !> @brief Initialises the id_r32_operator_pair object
  !> @param[in] operator_in  The operator that will be stored in the paired object
  !> @param[in] id           The integer ID to pair with the operator
  subroutine initialise( self, operator_in, id )

    implicit none

    class(id_r32_operator_pair_type), intent(inout) :: self
    type(operator_r32_type),          intent(in)    :: operator_in
    integer(kind=i_def),              intent(in)    :: id

    self%operator_ = operator_in%deep_copy()
    call self%set_id(id)

  end subroutine initialise

  !> @brief Get the operator corresponding to the paired object
  !> @param[in] self     The paired object
  !> @return             The operator
  function get_operator(self) result(operator_out)

    implicit none

    class(id_r32_operator_pair_type), target, intent(in) :: self
    type(operator_r32_type),                  pointer    :: operator_out

    operator_out => self%operator_

  end function get_operator

end module id_r32_operator_pair_mod
