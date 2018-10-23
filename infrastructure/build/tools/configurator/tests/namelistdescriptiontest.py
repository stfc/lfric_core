#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

import unittest
import StringIO
import random
import os
import tempfile

from subprocess import Popen

import configurator.namelistdescription as description

PICKER_EXE = 'rose_picker'


###############################################################################
class NamelistMetaTest(unittest.TestCase):

    def setUp(self):
        self.maxDiff = None


###############################################################################
    def tearDown(self):
        pass


###############################################################################
    def test_module_write_empty(self):
        output_file = StringIO.StringIO()

        uut = description.NamelistDescription('test')
        self.assertRaises(description.NamelistDescriptionException,
                          uut.write_module, output_file)


###############################################################################
    def test_module_write_one_of_each(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only : i_def, &
                            i_long, &
                            i_native, &
                            i_short, &
                            l_def, &
                            r_def, &
                            r_double, &
                            r_single, &
                            str_def, &
                            str_max_filename
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: enum_from_key, key_from_enum, &
            read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  integer(i_native), public, parameter :: test_enum_one = 135
  integer(i_native), public, parameter :: test_enum_three = 763
  integer(i_native), public, parameter :: test_enum_two = 847

  integer(i_def), public, protected :: dint
  logical(l_def), public, protected :: dlog
  real(r_def), public, protected :: dreal
  character(str_def), public, protected :: dstr
  integer(i_native), public, protected :: enum
  character(str_max_filename), public, protected :: fstr
  integer(i_long), public, protected :: lint
  real(r_double), public, protected :: lreal
  integer(i_short), public, protected :: sint
  real(r_single), public, protected :: sreal
  integer(i_def), public, protected :: vint
  real(r_def), public, protected :: vreal
  character(str_def), public, protected :: vstr

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: enum_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_native), parameter :: enum_value(3) &
          = [135_i_native, &
             763_i_native, &
             847_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function enum_from_key( key )

    use constants_mod, only: emdi
    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == emdi) then
      write( log_scratch_space, '(A)') &
          'Missing key for enum enumeration in test namelist.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_index = 1
    do
      if (trim(enum_key(key_index)) == trim(key)) then
        enum_from_key = enum_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(enum_key, 1)) then
          write( log_scratch_space, &
                 '("Key ''", A, "'' not recognised for test enum")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function enum_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_enum( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (enum_value(value_index) == value) then
        key_from_enum = enum_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(enum_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in test enum")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_enum

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_test_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank, enum )
  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, dummy_enum )

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_enum

    character(str_def) :: buffer_character_str_def(2)
    character(str_max_filename) :: buffer_character_str_max_filename(1)
    integer(i_def) :: buffer_integer_i_def(2)
    integer(i_long) :: buffer_integer_i_long(1)
    integer(i_native) :: buffer_integer_i_native(1)
    integer(i_short) :: buffer_integer_i_short(1)
    integer(i_native) :: buffer_logical_l_def(1)
    real(r_def) :: buffer_real_r_def(2)
    real(r_double) :: buffer_real_r_double(1)
    real(r_single) :: buffer_real_r_single(1)

    character(str_def) :: enum

    namelist /test/ dint, &
                    dlog, &
                    dreal, &
                    dstr, &
                    enum, &
                    fstr, &
                    lint, &
                    lreal, &
                    sint, &
                    sreal, &
                    vint, &
                    vreal, &
                    vstr

    integer(i_native) :: condition

    dint = imdi
    dlog = .false.
    dreal = rmdi
    dstr = cmdi
    enum = emdi
    fstr = cmdi
    lint = imdi
    lreal = rmdi
    sint = imdi
    sreal = rmdi
    vint = imdi
    vreal = rmdi
    vstr = cmdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_enum = enum_from_key( enum )

    end if

    buffer_integer_i_def(2) = dint
    buffer_logical_l_def(1) = merge( 1, 0, dlog )
    buffer_real_r_def(2) = dreal
    buffer_character_str_def(2) = dstr
    buffer_integer_i_native(1) = dummy_enum
    buffer_character_str_max_filename(1) = fstr
    buffer_integer_i_long(1) = lint
    buffer_real_r_double(1) = lreal
    buffer_integer_i_short(1) = sint
    buffer_real_r_single(1) = sreal
    buffer_integer_i_def(1) = vint
    buffer_real_r_def(1) = vreal
    buffer_character_str_def(1) = vstr

    call broadcast( buffer_character_str_def, 2*str_def, 0 )
    call broadcast( buffer_character_str_max_filename, 1*str_max_filename, 0 )
    call broadcast( buffer_integer_i_def, 2, 0 )
    call broadcast( buffer_integer_i_long, 1, 0 )
    call broadcast( buffer_integer_i_native, 1, 0 )
    call broadcast( buffer_integer_i_short, 1, 0 )
    call broadcast( buffer_logical_l_def, 1, 0 )
    call broadcast( buffer_real_r_def, 2, 0 )
    call broadcast( buffer_real_r_double, 1, 0 )
    call broadcast( buffer_real_r_single, 1, 0 )

    dint = buffer_integer_i_def(2)
    dlog = buffer_logical_l_def(1) /= 0
    dreal = buffer_real_r_def(2)
    dstr = buffer_character_str_def(2)
    dummy_enum = buffer_integer_i_native(1)
    fstr = buffer_character_str_max_filename(1)
    lint = buffer_integer_i_long(1)
    lreal = buffer_real_r_double(1)
    sint = buffer_integer_i_short(1)
    sreal = buffer_real_r_single(1)
    vint = buffer_integer_i_def(1)
    vreal = buffer_real_r_def(1)
    vstr = buffer_character_str_def(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    dint = imdi
    dlog = .false.
    dreal = real(rmdi,r_def)
    dstr = cmdi
    enum = int(imdi,i_native)
    fstr = cmdi
    lint = imdi
    lreal = real(rmdi,r_double)
    sint = imdi
    sreal = real(rmdi,r_single)
    vint = imdi
    vreal = real(rmdi,r_def)
    vstr = cmdi

    return
  end subroutine test_final


end module test_config_mod
        '''.strip()

        random.seed(1)
        uut = description.NamelistDescription('test')
        uut.add_value('vint', 'integer')
        uut.add_value('dint', 'integer', 'default')
        uut.add_value('sint', 'integer', 'short')
        uut.add_value('lint', 'integer', 'long')
        uut.add_value('dlog', 'logical', 'default')
        uut.add_value('vreal', 'real')
        uut.add_value('dreal', 'real', 'default')
        uut.add_value('sreal', 'real', 'single')
        uut.add_value('lreal', 'real', 'double')
        uut.add_string('vstr')
        uut.add_string('dstr', configure_string_length='default')
        uut.add_string('fstr', configure_string_length='filename')
        uut.add_enumeration('enum', enumerators=['one', 'two', 'three'])
        output_file = StringIO.StringIO()
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_module_write_growing(self):
        first_expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only : i_def, &
                            i_native
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  integer(i_def), public, protected :: foo

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_test_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank )
  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    integer(i_def) :: buffer_integer_i_def(1)

    namelist /test/ foo

    integer(i_native) :: condition

    foo = imdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_def(1) = foo

    call broadcast( buffer_integer_i_def, 1, 0 )

    foo = buffer_integer_i_def(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    foo = imdi

    return
  end subroutine test_final


end module test_config_mod
        '''.strip()

        second_expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the test namelist.
!>
module test_config_mod

  use constants_mod, only : i_def, &
                            i_native, &
                            r_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: read_test_namelist, postprocess_test_namelist, &
            test_is_loadable, test_is_loaded, test_final

  real(r_def), public, protected :: bar
  integer(i_def), public, protected :: foo

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_test_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank )
  end subroutine read_test_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    integer(i_def) :: buffer_integer_i_def(1)
    real(r_def) :: buffer_real_r_def(1)

    namelist /test/ bar, &
                    foo

    integer(i_native) :: condition

    bar = rmdi
    foo = imdi

    if (local_rank == 0) then

      read( file_unit, nml=test, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = bar
    buffer_integer_i_def(1) = foo

    call broadcast( buffer_integer_i_def, 1, 0 )
    call broadcast( buffer_real_r_def, 1, 0 )

    bar = buffer_real_r_def(1)
    foo = buffer_integer_i_def(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_test_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none


  end subroutine postprocess_test_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function test_is_loadable()

    implicit none

    logical :: test_is_loadable

    test_is_loadable = .not. namelist_loaded

  end function test_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function test_is_loaded()

    implicit none

    logical :: test_is_loaded

    test_is_loaded = namelist_loaded

  end function test_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine test_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    bar = real(rmdi,r_def)
    foo = imdi

    return
  end subroutine test_final


end module test_config_mod
        '''.strip()

        output_file = StringIO.StringIO()

        uut = description.NamelistDescription('test')
        uut.add_value('foo', 'integer')
        uut.write_module(output_file)

        self.assertMultiLineEqual(first_expected_source + '\n',
                                  output_file.getvalue())

        output_file = StringIO.StringIO()
        uut.add_value('bar', 'real', 'default')
        uut.write_module(output_file)

        self.assertMultiLineEqual(second_expected_source + '\n',
                                  output_file.getvalue())

    ##########################################################################
    def test_enumeration_only(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the enum namelist.
!>
module enum_config_mod

  use constants_mod, only : i_native, &
                            str_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: value_from_key, key_from_value, &
            read_enum_namelist, postprocess_enum_namelist, &
            enum_is_loadable, enum_is_loaded, enum_final

  integer(i_native), public, parameter :: enum_value_one = 135
  integer(i_native), public, parameter :: enum_value_three = 763
  integer(i_native), public, parameter :: enum_value_two = 847

  integer(i_native), public, protected :: value

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: value_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']

  integer(i_native), parameter :: value_value(3) &
          = [135_i_native, &
             763_i_native, &
             847_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function value_from_key( key )

    use constants_mod, only: emdi
    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == emdi) then
      write( log_scratch_space, '(A)') &
          'Missing key for value enumeration in enum namelist.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_index = 1
    do
      if (trim(value_key(key_index)) == trim(key)) then
        value_from_key = value_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(value_key, 1)) then
          write( log_scratch_space, &
                 '("Key ''", A, "'' not recognised for enum value")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function value_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_value( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (value_value(value_index) == value) then
        key_from_value = value_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(value_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in enum value")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_value

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_enum_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank, value )
  end subroutine read_enum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, dummy_value )

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_value

    integer(i_native) :: buffer_integer_i_native(1)

    character(str_def) :: value

    namelist /enum/ value

    integer(i_native) :: condition

    value = emdi

    if (local_rank == 0) then

      read( file_unit, nml=enum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_value = value_from_key( value )

    end if

    buffer_integer_i_native(1) = dummy_value

    call broadcast( buffer_integer_i_native, 1, 0 )

    dummy_value = buffer_integer_i_native(1)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_enum_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none


  end subroutine postprocess_enum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function enum_is_loadable()

    implicit none

    logical :: enum_is_loadable

    enum_is_loadable = .not. namelist_loaded

  end function enum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function enum_is_loaded()

    implicit none

    logical :: enum_is_loaded

    enum_is_loaded = namelist_loaded

  end function enum_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine enum_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    value = int(imdi,i_native)

    return
  end subroutine enum_final


end module enum_config_mod
        '''.strip()

        random.seed(1)
        uut = description.NamelistDescription('enum')
        uut.add_enumeration('value', enumerators=['one', 'two', 'three'])
        output_file = StringIO.StringIO()
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_more_than_one_enumeration(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the twoenum namelist.
!>
module twoenum_config_mod

  use constants_mod, only : i_native, &
                            str_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: first_from_key, key_from_first, &
            second_from_key, key_from_second, &
            read_twoenum_namelist, postprocess_twoenum_namelist, &
            twoenum_is_loadable, twoenum_is_loaded, twoenum_final

  integer(i_native), public, parameter :: twoenum_first_one = 135
  integer(i_native), public, parameter :: twoenum_first_three = 763
  integer(i_native), public, parameter :: twoenum_first_two = 847
  integer(i_native), public, parameter :: twoenum_second_ay = 256
  integer(i_native), public, parameter :: twoenum_second_bee = 495
  integer(i_native), public, parameter :: twoenum_second_see = 449

  integer(i_native), public, protected :: first
  integer(i_native), public, protected :: second

  logical :: namelist_loaded = .false.

  character(str_def), parameter :: first_key(3) &
          = [character(len=str_def) :: 'one', &
                                       'three', &
                                       'two']
  character(str_def), parameter :: second_key(3) &
          = [character(len=str_def) :: 'ay', &
                                       'bee', &
                                       'see']

  integer(i_native), parameter :: first_value(3) &
          = [135_i_native, &
             763_i_native, &
             847_i_native]
  integer(i_native), parameter :: second_value(3) &
          = [256_i_native, &
             495_i_native, &
             449_i_native]

contains

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function first_from_key( key )

    use constants_mod, only: emdi
    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == emdi) then
      write( log_scratch_space, '(A)') &
          'Missing key for first enumeration in twoenum namelist.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_index = 1
    do
      if (trim(first_key(key_index)) == trim(key)) then
        first_from_key = first_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(first_key, 1)) then
          write( log_scratch_space, &
                 '("Key ''", A, "'' not recognised for twoenum first")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function first_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_first( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (first_value(value_index) == value) then
        key_from_first = first_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(first_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in twoenum first")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_first

  !> Gets the enumeration value from the key string.
  !>
  !> An error is reported if the key is not actually a key.
  !>
  !> @param[in] key Enumeration key.
  !>
  integer(i_native) function second_from_key( key )

    use constants_mod, only: emdi
    implicit none

    character(*), intent(in) :: key

    integer(i_native) :: key_index

    if (key == emdi) then
      write( log_scratch_space, '(A)') &
          'Missing key for second enumeration in twoenum namelist.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    key_index = 1
    do
      if (trim(second_key(key_index)) == trim(key)) then
        second_from_key = second_value(key_index)
        return
      else
        key_index = key_index + 1
        if (key_index > ubound(second_key, 1)) then
          write( log_scratch_space, &
                 '("Key ''", A, "'' not recognised for twoenum second")' ) key
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function second_from_key

  !> Gets the enumeration key corresponding to a particular value.
  !>
  !> An error is reported if the value is not within range.
  !>
  !> @param[in] value Enumeration value.
  !>
  character(str_def) function key_from_second( value )

    implicit none

    integer(i_native), intent(in) :: value

    integer(i_native) :: value_index

    value_index = 1
    do
      if (second_value(value_index) == value) then
        key_from_second = second_key(value_index)
        return
      else
        value_index = value_index + 1
        if (value_index > ubound(second_key, 1)) then
          write( log_scratch_space, &
                 '("Value ", I0, " is not in twoenum second")' ) value
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
    end do

  end function key_from_second

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_twoenum_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank, first, second )
  end subroutine read_twoenum_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank, dummy_first, dummy_second )

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    integer(i_native), intent(out) :: dummy_first
    integer(i_native), intent(out) :: dummy_second

    integer(i_native) :: buffer_integer_i_native(2)

    character(str_def) :: first
    character(str_def) :: second

    namelist /twoenum/ first, &
                       second

    integer(i_native) :: condition

    first = emdi
    second = emdi

    if (local_rank == 0) then

      read( file_unit, nml=twoenum, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      dummy_first = first_from_key( first )
      dummy_second = second_from_key( second )

    end if

    buffer_integer_i_native(1) = dummy_first
    buffer_integer_i_native(2) = dummy_second

    call broadcast( buffer_integer_i_native, 2, 0 )

    dummy_first = buffer_integer_i_native(1)
    dummy_second = buffer_integer_i_native(2)



    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_twoenum_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none


  end subroutine postprocess_twoenum_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function twoenum_is_loadable()

    implicit none

    logical :: twoenum_is_loadable

    twoenum_is_loadable = .not. namelist_loaded

  end function twoenum_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function twoenum_is_loaded()

    implicit none

    logical :: twoenum_is_loaded

    twoenum_is_loaded = namelist_loaded

  end function twoenum_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine twoenum_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    first = int(imdi,i_native)
    second = int(imdi,i_native)

    return
  end subroutine twoenum_final


end module twoenum_config_mod
                       '''.strip()

        random.seed(1)
        uut = description.NamelistDescription('twoenum')
        uut.add_enumeration('first', enumerators=['one', 'two', 'three'])
        uut.add_enumeration('second', enumerators=['ay', 'bee', 'see'])
        output_file = StringIO.StringIO()
        uut.write_module(output_file)

        self.assertMultiLineEqual(output_file.getvalue(),
                                  expected_source + '\n')

    ###########################################################################
    def test_module_write_computed(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the teapot namelist.
!>
module teapot_config_mod

  use constants_mod, only : i_native, &
                            r_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: read_teapot_namelist, postprocess_teapot_namelist, &
            teapot_is_loadable, teapot_is_loaded, teapot_final

  real(r_def), public, protected :: bar
  real(r_def), public, protected :: foo

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_teapot_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank )
  end subroutine read_teapot_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    real(r_def) :: buffer_real_r_def(1)

    namelist /teapot/ foo

    integer(i_native) :: condition

    bar = rmdi
    foo = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=teapot, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = foo

    call broadcast( buffer_real_r_def, 1, 0 )

    foo = buffer_real_r_def(1)


    bar = foo ** 2

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_teapot_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi

    implicit none


  end subroutine postprocess_teapot_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function teapot_is_loadable()

    implicit none

    logical :: teapot_is_loadable

    teapot_is_loadable = .not. namelist_loaded

  end function teapot_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function teapot_is_loaded()

    implicit none

    logical :: teapot_is_loaded

    teapot_is_loaded = namelist_loaded

  end function teapot_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine teapot_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    bar = real(rmdi,r_def)
    foo = real(rmdi,r_def)

    return
  end subroutine teapot_final


end module teapot_config_mod
        '''.strip()

        output_file = StringIO.StringIO()

        uut = description.NamelistDescription('teapot')
        uut.add_value('foo', 'real', 'default')
        uut.add_computed('bar', 'real', 'default', calculation=['foo ** 2'])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_module_write_constant(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the cheese namelist.
!>
module cheese_config_mod

  use constants_mod, only : i_native, &
                            r_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: read_cheese_namelist, postprocess_cheese_namelist, &
            cheese_is_loadable, cheese_is_loaded, cheese_final

  real(r_def), public, protected :: fred
  real(r_def), public, protected :: wilma

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_cheese_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank )
  end subroutine read_cheese_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only : cmdi, emdi, FUDGE, imdi, rmdi

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    real(r_def) :: buffer_real_r_def(1)

    namelist /cheese/ fred

    integer(i_native) :: condition

    fred = rmdi
    wilma = rmdi

    if (local_rank == 0) then

      read( file_unit, nml=cheese, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_real_r_def(1) = fred

    call broadcast( buffer_real_r_def, 1, 0 )

    fred = buffer_real_r_def(1)


    wilma = fred * FUDGE

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_cheese_namelist()

    use constants_mod, only : cmdi, emdi, FUDGE, imdi, rmdi

    implicit none


  end subroutine postprocess_cheese_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function cheese_is_loadable()

    implicit none

    logical :: cheese_is_loadable

    cheese_is_loadable = .not. namelist_loaded

  end function cheese_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function cheese_is_loaded()

    implicit none

    logical :: cheese_is_loaded

    cheese_is_loaded = namelist_loaded

  end function cheese_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine cheese_final()

    use constants_mod, only: cmdi, emdi, FUDGE, imdi, rmdi

    implicit none

    fred = real(rmdi,r_def)
    wilma = real(rmdi,r_def)

    return
  end subroutine cheese_final


end module cheese_config_mod
        '''.strip()

        output_file = StringIO.StringIO()

        uut = description.NamelistDescription('cheese')
        uut.add_usage('FUDGE', module='constants_mod')
        uut.add_value('fred', 'real', 'default')
        uut.add_computed('wilma', 'real', 'default',
                         calculation=['fred * FUDGE'])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_module_write_array(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> Manages the aerial namelist.
!>
module aerial_config_mod

  use constants_mod, only : i_def, &
                            i_native, &
                            r_def, &
                            str_def
  use log_mod,       only : log_event, log_scratch_space, LOG_LEVEL_ERROR
  use mpi_mod,       only : broadcast
  use mpi,           only : MPI_SUCCESS

  implicit none

  private
  public :: read_aerial_namelist, postprocess_aerial_namelist, &
            aerial_is_loadable, aerial_is_loaded, aerial_final

  integer(i_native), parameter, public :: max_array_size = 100

  character(str_def), public, protected :: absolute(5)
  integer(i_def), public, protected, allocatable :: inlist(:)
  integer(i_native), public, protected :: lsize
  real(r_def), public, protected, allocatable :: outlist(:)
  integer(i_def), public, protected, allocatable :: unknown(:)

  logical :: namelist_loaded = .false.

contains

  !> Populates this module from a namelist file.
  !>
  !> An error is reported if the namelist could not be read.
  !>
  !> @param [in] file_unit Unit number of the file to read from.
  !> @param [in] local_rank Rank of current process.
  !>
  subroutine read_aerial_namelist( file_unit, local_rank )
    implicit none
    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank
    call read_namelist( file_unit, local_rank )
  end subroutine read_aerial_namelist

  ! Reads the namelist file.
  !
  subroutine read_namelist( file_unit, local_rank )

    use constants_mod, only : cmdi, emdi, imdi, rmdi
    use wibble_mod, only : esize

    implicit none

    integer(i_native), intent(in) :: file_unit
    integer(i_native), intent(in) :: local_rank

    integer(i_native) :: buffer_integer_i_native(1)

    namelist /aerial/ absolute, &
                      inlist, &
                      lsize, &
                      outlist, &
                      unknown

    integer(i_native) :: condition

    allocate( inlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( outlist(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    allocate( unknown(max_array_size), stat=condition )
    if (condition /= 0) then
      write( log_scratch_space, '(A)' ) &
            'Unable to allocate temporary array for "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    absolute = cmdi
    inlist = imdi
    lsize = imdi
    outlist = rmdi
    unknown = imdi

    if (local_rank == 0) then

      read( file_unit, nml=aerial, iostat=condition, iomsg=log_scratch_space )
      if (condition /= 0) then
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

    end if

    buffer_integer_i_native(1) = lsize

    call broadcast( buffer_integer_i_native, 1, 0 )

    lsize = buffer_integer_i_native(1)



    call broadcast( absolute, size(absolute, 1)*str_def, 0 )
    call broadcast( inlist, size(inlist, 1), 0 )
    call broadcast( outlist, size(outlist, 1), 0 )
    call broadcast( unknown, size(unknown, 1), 0 )

    namelist_loaded = .true.

  end subroutine read_namelist

  !> Performs any processing to be done once all namelists are loaded
  !>
  subroutine postprocess_aerial_namelist()

    use constants_mod, only : cmdi, emdi, imdi, rmdi
    use wibble_mod, only : esize

    implicit none

    integer(i_native) :: condition
    integer(i_native) :: index
    integer(i_native) :: size
    integer(i_def), allocatable :: new_inlist(:)
    real(r_def), allocatable :: new_outlist(:)
    integer(i_def), allocatable :: new_unknown(:)

    size = lsize
    allocate( new_inlist(size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "inlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_inlist(:size) = inlist(:size)
    call move_alloc( new_inlist, inlist )
    if (allocated(new_inlist)) deallocate( new_inlist)
    size = esize
    allocate( new_outlist(size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "outlist"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_outlist(:size) = outlist(:size)
    call move_alloc( new_outlist, outlist )
    if (allocated(new_outlist)) deallocate( new_outlist)
    do index=ubound(unknown, 1), 1, -1
      if (unknown(index) /= imdi) exit
    end do
    size = index
    allocate( new_unknown(size), stat=condition )
    if (condition /= 0) then
      write(log_scratch_space, '(A)') 'Unable to allocate "unknown"'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    new_unknown(:size) = unknown(:size)
    call move_alloc( new_unknown, unknown )
    if (allocated(new_unknown)) deallocate( new_unknown)

  end subroutine postprocess_aerial_namelist

  !> Can this namelist be loaded?
  !>
  !> @return True if it is possible to load the namelist.
  !>
  function aerial_is_loadable()

    implicit none

    logical :: aerial_is_loadable

    aerial_is_loadable = .not. namelist_loaded

  end function aerial_is_loadable

  !> Has this namelist been loaded?
  !>
  !> @return True if the namelist has been loaded.
  !>
  function aerial_is_loaded()

    implicit none

    logical :: aerial_is_loaded

    aerial_is_loaded = namelist_loaded

  end function aerial_is_loaded

  !> Clear out any allocated memory
  !>
  subroutine aerial_final()

    use constants_mod, only: cmdi, emdi, imdi, rmdi

    implicit none

    absolute = cmdi
    lsize = imdi

    if ( allocated(inlist) ) deallocate(inlist)
    if ( allocated(outlist) ) deallocate(outlist)
    if ( allocated(unknown) ) deallocate(unknown)

    return
  end subroutine aerial_final


end module aerial_config_mod
        '''.strip()

        output_file = StringIO.StringIO()

        uut = description.NamelistDescription('aerial')
        uut.add_usage('esize', module='wibble_mod')
        uut.add_value('lsize', 'integer', 'native')
        uut.add_string('absolute', bounds=['5'])
        uut.add_value('inlist', 'integer', bounds=['lsize'])
        uut.add_value('outlist', 'real', bounds=['esize'])
        uut.add_value('unknown', 'integer', bounds=[':'])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())


##############################################################################
class NamelistConfigDescriptionTest(unittest.TestCase):
    ##########################################################################
    def setUp(self):
        self.nml_config_file = None

    ##########################################################################
    def tearDown(self):
        if self.nml_config_file:
            os.remove(self.nml_config_file)

        if os.path.isfile('config_namelists.txt'):
            os.remove('config_namelists.txt')

    ##########################################################################
    @staticmethod
    def _compile_dictionary(namelists):

        dictionary = {}

        for namelist in namelists:

            parameters = {}
            for parameter in namelist.get_parameters():

                parameters[parameter.name] = \
                    [parameter.fortran_type.intrinsic_type,
                     parameter.fortran_type.kind]
                if parameter.get_configure_type() == 'enumeration':
                    parameters[parameter.name].extend(
                        parameter.mapping.keys())
                elif parameter.get_configure_type() == 'computed':
                    parameters[parameter.name].append(
                        parameter.computation)
                elif parameter.get_configure_type() == 'array':
                    if isinstance(parameter.bounds, list):
                        bounds = parameter.bounds[0]
                    else:
                        bounds = parameter.bounds
                    parameters[parameter.name].append('(' + bounds + ')')

            dictionary[namelist.get_namelist_name()] = parameters

        return dictionary

    ##########################################################################
    def test_parser_good_file(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:fred]

[namelist:fred=first_thing]
type=character

[namelist:fred=second]
type=integer

[namelist:fred=filename]
type=character
!string_length=filename

[namelist:fred=choices]
!enumeration=true
values=foo, bar, baz, qux
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(self.nml_config_file))

        self.assertEqual({'fred':
                          {'choices':     ['integer', 'i_native',
                                           'foo', 'bar', 'baz', 'qux'],
                           'filename':    ['character', 'str_max_filename'],
                           'first_thing': ['character', 'str_def'],
                           'second':      ['integer', 'i_def']}},
                         result)

    ##########################################################################
    def test_only_enumeration(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:barney]

[namelist:barney=stuff]
!enumeration=true
values=one, two, three
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(self.nml_config_file))

        self.assertEqual({'barney': {'stuff': ['integer', 'i_native',
                                               'one', 'two', 'three']}},
                         result)

    ##########################################################################
    def test_non_enumeration_no_type(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:barney]

[namelist:barney=stuff]
values=one, two, three
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        self.assertRaises(description.NamelistDescriptionException,
                          uut.process_config, self.nml_config_file)

    ##########################################################################
    def test_no_member_type(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:santa]

[namelist:santa=elf]
length=:
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        self.assertRaises(description.NamelistDescriptionException,
                          uut.process_config, self.nml_config_file)

    ##########################################################################
    def test_computed_fields(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:teapot]

[namelist:teapot=foo]
type=real

[!namelist:teapot=bar]
type=real
expression=namelist:teapot=foo ** 2

[!namelist:teapot=baz]
type=real
expression=source:constants_mod=PI * foo
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(self.nml_config_file))

        self.assertEqual({'teapot': {'foo': ['real', 'r_def'],
                                     'bar': ['real', 'r_def', 'foo ** 2'],
                                     'baz': ['real', 'r_def', 'PI * foo']}},
                         result)


##########################################################################
    def test_constant_in_computed(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:cheese]

[namelist:cheese=fred]
type=real

[!namelist:cheese=wilma]
type=real
expression=fred * source:constants_mod=FUDGE
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(self.nml_config_file))

        self.assertEqual(
            {'cheese': {'fred':  ['real', 'r_def'],
                        'wilma': ['real', 'r_def', 'fred * FUDGE']}},
            result)


##########################################################################
    def test_array_fields(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:aerial]

[namelist:aerial=fred]
type=real

[namelist:aerial=wilma]
type=real
length=:
!bounds=source:constants_mod=FUDGE

[namelist:aerial=betty]
type=logical
length=:
!bounds=fred

[namelist:aerial=dino]
type=integer
length=:
!bounds=namelist:sugar=TABLET

[namelist:aerial=bambam]
type=integer
length=:
''')
        input_file.seek(0)

        picker_command = "{} {}".format(PICKER_EXE, input_file.name)
        out = Popen(picker_command, shell=True)
        out.wait()

        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))
        input_file.close()

        uut = description.NamelistConfigDescription()
        result = self._compile_dictionary(
            uut.process_config(self.nml_config_file))

        self.assertEqual(
            {'aerial': {'bambam': ['integer', 'i_def', '(:)'],
                        'betty':  ['logical', 'l_def', '(fred)'],
                        'fred':   ['real', 'r_def'],
                        'wilma':  ['real', 'r_def', '(FUDGE)'],
                        'dino':   ['integer', 'i_def', '(TABLET)']}},
            result)


##############################################################################
if __name__ == '__main__':
    unittest.main()
