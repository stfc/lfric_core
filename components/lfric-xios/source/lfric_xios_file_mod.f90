!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module containing a file type for XIOS interface
!>
module lfric_xios_file_mod

  use constants_mod,         only: i_def, str_def, str_max_filename
  use linked_list_data_mod,  only: linked_list_data_type

implicit none

private

!> @brief Container for file properties need by XIOS
!>
type, extends(linked_list_data_type), public :: xios_file_type
  private

  !> Unique identifier for XIOS file handle
  character(str_def) :: xios_id
  !> Path to file
  character(str_max_filename) :: path
  !> File output frequency
  integer(i_def) :: output_freq
  !> XIOS ID of associated field group
  character(str_def) :: field_group = "unset"

contains
  procedure, public :: init_xios_file
  procedure, public :: get_xios_id
  procedure, public :: get_path
  procedure, public :: get_output_freq
  procedure, public :: get_field_group

end type xios_file_type

contains
!-------------------------------------------------------------------------------
! XIOS file type procedures
!-------------------------------------------------------------------------------

!>  @brief  Initialises a file object for use with XIOS
!>
!>  @param[in]           input_xios_id   The XIOS name for the file
!>  @param[in,optional]  path            The path to the file
!>  @param[in,optional]  freq            The file output frequency
!>  @param[in,optional]  field_group_id  The associated field group id
!>
subroutine init_xios_file(self, input_xios_id, path, freq, field_group_id)

  implicit none

  class(xios_file_type),      intent(inout) :: self
  character(len=*),           intent(in)    :: input_xios_id
  character(len=*), optional, intent(in)    :: path
  integer(i_def),   optional, intent(in)    :: freq
  character(len=*), optional, intent(in)    :: field_group_id

  self%xios_id = input_xios_id

  ! If path is not present then default to use xios_id as file name in working
  ! directory
  if (present(path)) then
    self%path = path
  else
    self%path = input_xios_id
  end if

  if (present(freq)) then
    self%output_freq = freq
  else
    ! -999 used as flag for frequency not to be set
    self%output_freq = -999
  end if

  if (present(field_group_id)) self%field_group = field_group_id

end subroutine init_xios_file

!> Getter for XIOS file ID
!> @param[out]  output_xios_id  The XIOS ID for the file
!>
function get_xios_id(self) result(output_xios_id)

  implicit none

  class(xios_file_type),   intent(inout) :: self

  character(str_def) :: output_xios_id

  output_xios_id = self%xios_id

end function get_xios_id

!> Getter for file path
!> @param[out] output_path The path to the file
!>
function get_path(self) result(output_path)

  implicit none

  class(xios_file_type),            intent(inout) :: self

  character(str_max_filename) :: output_path

  output_path = self%path

end function get_path

!> Getter for file output frequency
!> @param[out]  file_freq  The file output frequency
!>
function get_output_freq(self) result(file_freq)

  implicit none

  class(xios_file_type), intent(inout) :: self

  integer(i_def) :: file_freq

  file_freq = self%output_freq

end function get_output_freq

!> Getter for associated field group ID
!> @param[out]  field_group_id  The associated field group id
!>
function get_field_group(self) result(field_group_id)

  implicit none

  class(xios_file_type), intent(inout) :: self

  character(str_def) :: field_group_id

  field_group_id = self%field_group

end function get_field_group

end module lfric_xios_file_mod
