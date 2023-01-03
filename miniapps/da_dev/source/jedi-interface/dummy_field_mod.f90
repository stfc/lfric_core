!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a mock JEDI field class.
!>
!> @details This module holds a mock JEDI field class that included to provide
!>          a simple field container for external field data. In JEDI we are
!>          using Atlas fields. The data is stored continuous vertically and
!>          unstructured horizontally. The data is a 2D array where the first
!>          index stores columns and the second the horizontal.
!>
module dummy_field_mod

  use, intrinsic :: iso_fortran_env, only : real64
  use fs_continuity_mod,             only : W3, Wtheta
  use constants_mod,                 only : i_def, str_def
  use log_mod,                       only : log_event, log_scratch_space, &
                                            LOG_LEVEL_ERROR

  implicit none

  private

type, public :: dummy_field_type
  private

  !> The 64-bit floating point values of the field
  real( kind=real64 ), allocatable   :: data(:,:)
  !> the name of the field
  character( len=str_def )           :: field_name
  !> number of vertical points in the external data
  integer( kind=i_def )              :: n_vertical
  !> number of horizontal points in the external data
  integer( kind=i_def )              :: n_horizontal

contains

  !> Field initialiser.
  procedure, public :: initialise => dummy_field_initialiser

  !> Get a pointer to the dummy field
  procedure, public :: get_data

  !> Get the external field name
  procedure, public :: get_field_name

  !> Finalizer
  final             :: dummy_field_destructor

end type dummy_field_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> dummy_field constructor
!> @param [in] n_vertical the number of vertical points
!> @param [in] n_horizontal the number of horizontal points
!> @param [in] fs function space enumerator
!> @param [in] is_2d logical defining if the field is 2d
!> @param [in] field_name name of the field
subroutine dummy_field_initialiser( self, n_vertical, n_horizontal, fs, is_2d, field_name )

  implicit none

  class( dummy_field_type ), intent(inout)    :: self
  integer( kind=i_def ), intent(in)           :: n_vertical
  integer( kind=i_def ), intent(in)           :: n_horizontal
  integer( kind=i_def ), intent(in)           :: fs
  logical, intent(in)                         :: is_2d
  character( len=* ), optional, intent(in)    :: field_name

  ! Local
  integer( kind=i_def )                       :: n_vertical_local

  ! setup
  if (is_2d) then
    n_vertical_local = 1
  elseif (fs==W3) then
    n_vertical_local = n_vertical
  elseif (fs==Wtheta) then
    n_vertical_local = n_vertical + 1
  else
    log_scratch_space = 'The requested field type is not supported.'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  endif

  allocate(self%data(n_vertical_local,n_horizontal))
  self % field_name = field_name
  self % n_vertical = n_vertical_local
  self % n_horizontal = n_horizontal

end subroutine dummy_field_initialiser

!> Get a pointer to the data
!> @param [in] data_ptr a 2d real array to the field data
subroutine get_data(self, data_ptr)

  implicit none

  class( dummy_field_type ), target, intent(inout) :: self
  real(real64), pointer, intent(inout)             :: data_ptr(:,:)

  data_ptr => self % data

end subroutine get_data

!> Returns the external name of the field
!> @param [out] field_name the name of the field
function get_field_name( self ) result( field_name )

  implicit none

  class( dummy_field_type ), intent(in) :: self
  character( len=str_def )              :: field_name

  field_name = self%field_name

end function get_field_name

!> dummy_field finalizer
subroutine dummy_field_destructor(self)!

  implicit none

  type(dummy_field_type), intent(inout)    :: self

  if ( allocated(self%data) ) deallocate(self%data)
  self % field_name = ""

end subroutine dummy_field_destructor

end module dummy_field_mod
