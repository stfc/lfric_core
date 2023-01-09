!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Base class for paired objects, which describe fields,
!!        operators or other entities that are paired with a unique identifier.
module id_abstract_pair_mod

  use constants_mod,        only: i_def
  use linked_list_data_mod, only: linked_list_data_type

  implicit none

  private

  ! ========================================================================== !
  ! Base ID-Something Pair Object
  ! ========================================================================== !

  ! Public types
  type, extends(linked_list_data_type), public, abstract :: id_abstract_pair_type

    private

    integer(i_def), allocatable :: gnu_dummy(:)
    integer(i_def) :: id_

  contains

    procedure, public :: get_id
    procedure, public :: set_id

  end type id_abstract_pair_type

contains

  !> @brief Get the integer ID corresponding to the paired object
  !> @param[in] self  The paired object
  !> @return          The id
  function get_id(self) result(id)

    implicit none

    class(id_abstract_pair_type), intent(in) :: self
    integer(kind=i_def)                      :: id

    id = self%id_

  end function get_id

  !> @brief Sets the integer ID for the paired object
  !> @param[in] self  The paired object
  !> @param[in] id    The integer ID
  subroutine set_id(self, id)

    implicit none

    class(id_abstract_pair_type), intent(inout) :: self
    integer(kind=i_def),          intent(in)    :: id

    self%id_ = id

  end subroutine set_id

end module id_abstract_pair_mod
