!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a mock JEDI Geometry class.
!>
!> @details This module holds a mock JEDI Geometry class that includes only the
!>          functionality that is required to support the mock JEDI interface.
!>
module jedi_geometry_mod

  use, intrinsic :: iso_fortran_env, only : real64

  use log_mod,                       only : log_event,          &
                                            log_scratch_space,  &
                                            LOG_LEVEL_ERROR
  use constants_mod,                 only : i_def, str_def, l_def

  implicit none

  private

type, public :: jedi_geometry_type
  private

  !> The data map between external field data and LFRic fields
  integer( kind = i_def ), allocatable :: horizontal_map(:)
  !> the LFRic field dimensions
  integer( kind = i_def )              :: n_vertical
  integer( kind = i_def )              :: n_horizontal

contains

  !> Field initialiser.
  procedure, public :: initialise => jedi_geometry_initialiser

  !> getters
  procedure, public :: get_n_horizontal
  procedure, public :: get_n_vertical
  procedure, public :: get_horizontal_map

  !> Finalizer
  final             :: jedi_geometry_destructor

end type jedi_geometry_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> jedi_geometry constructor
subroutine jedi_geometry_initialiser( self )

  use da_dev_driver_mod, only : mesh

  implicit none

  class( jedi_geometry_type ), intent(inout)   :: self

  ! Local
  integer :: i_horizontal

  ! These will be provided by calls that James is working on
  self%n_vertical = mesh%get_nlayers()
  self%n_horizontal = mesh%get_last_edge_cell()

  ! Create horizontal_map
  allocate(self%horizontal_map(self%n_horizontal))

  do i_horizontal=1,self%n_horizontal
    self%horizontal_map(i_horizontal) = i_horizontal
  enddo

end subroutine jedi_geometry_initialiser

!> Get the number of horizontal points
function get_n_horizontal(self) result(n_horizontal)

  implicit none

  class( jedi_geometry_type ), intent(in) :: self
  integer(kind = i_def) :: n_horizontal

  n_horizontal = self%n_horizontal

end function get_n_horizontal

!> Get the number of vertical points
function get_n_vertical(self) result(n_vertical)

  implicit none

  class( jedi_geometry_type ), intent(in) :: self
  integer(kind = i_def)              :: n_vertical

  n_vertical = self%n_vertical

end function get_n_vertical

!> Get a pointer to the horizontal map
subroutine get_horizontal_map(self, horizontal_map)

  implicit none

  class( jedi_geometry_type ), target, intent(in) :: self
  integer(i_def), pointer, intent(inout)     :: horizontal_map(:)

  horizontal_map => self % horizontal_map

end subroutine get_horizontal_map

!> jedi_geometry finalizer
subroutine jedi_geometry_destructor(self)!

  implicit none

  type(jedi_geometry_type), intent(inout)    :: self

  if ( allocated(self % horizontal_map) ) deallocate(self % horizontal_map)

end subroutine jedi_geometry_destructor

end module jedi_geometry_mod
