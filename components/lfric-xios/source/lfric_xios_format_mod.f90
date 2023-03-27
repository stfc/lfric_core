!-------------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Field formatting for XIOS

module lfric_xios_format_mod

  use, intrinsic :: iso_fortran_env,  only : real32, real64
  use constants_mod,                  only: i_def
  use lfric_xios_constants_mod,       only: dp_xios, xios_max_int
  use field_parent_mod,               only: field_parent_proxy_type
  use field_r32_mod,                  only: field_r32_proxy_type
  use field_r64_mod,                  only: field_r64_proxy_type
  use integer_field_mod,              only: integer_field_proxy_type
  use log_mod,                        only:                                   &
                                        log_event,                            &
                                        log_level_warning,                    &
                                        log_level_error

  implicit none

  public :: format_field, inverse_format_field

  private
contains

!> @brief Formats LFRic field data array in accordance with XIOS conventions.
!> @details The LFRic field data form a matrix with the columns representing
!> vertical slices. To put it differently, denoting the matrix by a_ki,
!> the column index k selects the vertical level, whereas the row index i
!> selects the dof at ground level.
!> Since Fortran lays out matrices in column-major form ("contiguous
!> columns"), this makes iteration in the vertical direction
!> ("k-first") fast.
!> XIOS, on the other hand, expects the matrix in the form a_ik with the
!> horizontal index i preceding the vertical index k. Therefore, the
!> matrix has to be transposed before passing it to XIOS, and that's
!> what's happening here.
!> To understand the code, note that the data are passed in as a
!> one-dimensional array with the data laid out by columns, as already
!> observed.
!> @pre The total dimension must be an integer multiple of the number of rows.
!> @param[out] xios_data         Data array to be sent to XIOS
!> @param[in]  field_name        Field name for error reporting
!> @param[in]  fpxy              Field proxy of LFRic field to be formatted
!> @param[in]  m                 Number of rows in input matrix
!> @param[in]  n                 Number of columns in input matrix
!>
subroutine format_field(xios_data, field_name, fpxy, m, n)
  implicit none

  real(dp_xios),                  intent(in out):: xios_data(:)
  character(len=*),               intent(in)    :: field_name
  class(field_parent_proxy_type), intent(in)    :: fpxy
  integer(i_def),                 intent(in)    :: m
  integer(i_def),                 intent(in)    :: n

  integer(i_def) :: mn
  integer(i_def) :: i

  mn = m*n

  ! sanity check
  if (mn /= size(xios_data)) then
    call log_event('assertion failed for field ' // field_name                &
      // ': mn must equal size(xios_data)',  log_level_error)
  end if

  select type(fpxy)
    type is (field_r32_proxy_type)
      ! this elegant call can be up to 4 times slower than the loop below:
      ! xios_data = reshape(transpose(reshape(fpxy%data(:), (/m, n/))), (/mn/))
      do i = 1, m
        xios_data(1+(i-1)*n:i*n) = fpxy%data(i:mn:m)
      end do
    type is (field_r64_proxy_type)
      ! xios_data = reshape(transpose(reshape(fpxy%data(:), (/m, n/))), (/mn/))
      do i = 1, m
        xios_data(1+(i-1)*n:i*n) = fpxy%data(i:mn:m)
      end do

    type is (integer_field_proxy_type)
      if ( any( abs(fpxy%data(1:mn)) > xios_max_int) ) then
        call log_event( 'data for integer field ' //                          &
          trim(adjustl(field_name)) //                                        &
          ' contains values too large for 16-bit precision', log_level_warning)
      end if
      ! xios_data = reshape(transpose(reshape(fpxy%data(:), (/m, n/))), (/mn/))
      do i = 1, m
        xios_data(1+(i-1)*n:i*n) = fpxy%data(i:mn:m)
      end do

    class default
      call log_event( "invalid type for input field proxy", log_level_error )
  end select

end subroutine format_field

!> @brief Inverse of format_field above.
!> @param[out] xios_data         Data array received from XIOS
!> @param[in]  field_name        Field name for error reporting
!> @param[in]  fpxy              Field proxy of LFRic field to be read
!> @param[in]  m                 Number of rows in field data matrix
!> @param[in]  n                 Number of columns in field data matrix
!>
subroutine inverse_format_field(xios_data, field_name, fpxy, m, n)
  implicit none

  real(dp_xios),                  intent(in)     :: xios_data(:)
  character(len=*),               intent(in)     :: field_name
  class(field_parent_proxy_type), intent(in out) :: fpxy
  integer(i_def),                 intent(in)     :: m
  integer(i_def),                 intent(in)     :: n

  integer(i_def) :: mn
  integer(i_def) :: i

  mn = m*n

  ! sanity check
  if (mn /= size(xios_data)) then
    call log_event('assertion failed for field ' // field_name                &
      // ': mn must equal size(xios_data)',  log_level_error)
  end if

  ! Reshape the data to what we require for the LFRic field.
  ! Note the conversion from dp_xios to real32, real64 or i_def.
  select type(fpxy)
  type is (field_r32_proxy_type)
    do i = 1, m
      fpxy%data(i:mn:m) = real(xios_data(1+(i-1)*n:i*n), real32)
    end do

  type is (field_r64_proxy_type)
    do i = 1, m
      fpxy%data(i:mn:m) = real(xios_data(1+(i-1)*n:i*n), real64)
    end do

  type is (integer_field_proxy_type)
    do i = 1, m
      fpxy%data(i:mn:m) = int(xios_data(1+(i-1)*n:i*n), i_def)
    end do

  class default
    call log_event( "Invalid type for input field proxy", log_level_error)
  end select

end subroutine inverse_format_field

end module lfric_xios_format_mod
