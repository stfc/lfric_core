#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

import unittest
import StringIO

import configurator.namelistdescription as namelist
import configurator.namelistfeigner as feigner


###############################################################################
class FeignerTest(unittest.TestCase):
    def setUp(self):
        self.maxDiff = None

    ###########################################################################
    def test_empty(self):

        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module empty_mod

  use constants_mod, only : i_native
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

end module empty_mod
        '''.strip()

        output_file = StringIO.StringIO()
        uut = feigner.NamelistFeigner('empty_mod')
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_simple(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module simple_mod

  use constants_mod, only : i_def, i_native, l_def, r_double, str_def
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private
  public :: feign_simple_config

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_simple_config( foo, &
                                  bar, &
                                  baz, &
                                  fred )

    use simple_config_mod, only : read_simple_namelist

    implicit none

    integer(i_def), intent(in) :: foo
    real(r_double), intent(in) :: bar
    character(*), intent(in) :: baz
    logical(l_def), intent(in) :: fred

    integer(i_native) :: condition

    if (local_rank == -1) then
      local_rank = get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_simple_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&simple")' )
    write( temporary_unit, '("foo = ", I0)' ) foo
    write( temporary_unit, '("bar = ", E14.7)' ) bar
    write( temporary_unit, '("baz = ''", A, "''")' ) baz
    write( temporary_unit, '("fred = ", L2)' ) fred
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_simple_namelist( temporary_unit, local_rank )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop &
        'feign_simple_config: '// &
        'Unable to close temporary file'

  end subroutine feign_simple_config

end module simple_mod
        '''.strip()

        simple = namelist.NamelistDescription('simple')
        simple.add_value('foo', 'integer', 'default')
        simple.add_value('bar', 'real', 'double')
        simple.add_string('baz')
        simple.add_usage('qux', 'constants_mod')
        simple.add_value('fred', 'logical')

        output_file = StringIO.StringIO()
        uut = feigner.NamelistFeigner('simple_mod')
        uut.add_namelist([simple])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_enumeration(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module enumeration_mod

  use constants_mod, only : i_native
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private
  public :: feign_enum_config

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_enum_config( thing )

    use enum_config_mod, only : read_enum_namelist, &
                                key_from_thing, &
                                thing_from_key

    implicit none

    integer(i_native), intent(in) :: thing

    integer(i_native) :: condition

    if (local_rank == -1) then
      local_rank = get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_enum_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&enum")' )
    write( temporary_unit, '("thing = ''", A, "''")' ) key_from_thing( thing )
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_enum_namelist( temporary_unit, local_rank )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop &
        'feign_enum_config: '// &
        'Unable to close temporary file'

  end subroutine feign_enum_config

end module enumeration_mod
        '''.strip()

        enumable = namelist.NamelistDescription('enum')
        enumable.add_enumeration('thing', enumerators=['one', 'two'])

        output_file = StringIO.StringIO()
        uut = feigner.NamelistFeigner('enumeration_mod')
        uut.add_namelist([enumable])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_computed(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module computed_mod

  use constants_mod, only : i_def, i_native
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private
  public :: feign_computed_config

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_computed_config( teapot, &
                                    cheese )

    use computed_config_mod, only : read_computed_namelist

    implicit none

    integer(i_def), intent(in) :: teapot
    integer(i_def), intent(in) :: cheese

    integer(i_native) :: condition

    if (local_rank == -1) then
      local_rank = get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_computed_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&computed")' )
    write( temporary_unit, '("teapot = ", I0)' ) teapot
    write( temporary_unit, '("cheese = ", I0)' ) cheese
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_computed_namelist( temporary_unit, local_rank )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop &
        'feign_computed_config: '// &
        'Unable to close temporary file'

  end subroutine feign_computed_config

end module computed_mod
        '''.strip()

        simple = namelist.NamelistDescription('computed')
        simple.add_value('teapot', 'integer', 'default')
        simple.add_value('cheese', 'integer', 'default')
        simple.add_computed('biscuits', 'integer', 'default',
                            calculation=['teapot + cheese'])

        output_file = StringIO.StringIO()
        uut = feigner.NamelistFeigner('computed_mod')
        uut.add_namelist([simple])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_everything(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module everything_mod

  use constants_mod, only : i_native, l_def, r_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private
  public :: feign_everything_config

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_everything_config( cake, &
                                      teapot, &
                                      cheese, &
                                      fish, &
                                      tail )

    use everything_config_mod, only : read_everything_namelist, &
                                      key_from_teapot, &
                                      teapot_from_key

    implicit none

    character(*), intent(in) :: cake
    integer(i_native), intent(in) :: teapot
    logical(l_def), intent(in) :: cheese
    real(r_def), intent(in) :: fish
    integer(i_native), intent(in) :: tail

    integer(i_native) :: condition

    if (local_rank == -1) then
      local_rank = get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_everything_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&everything")' )
    write( temporary_unit, '("cake = ''", A, "''")' ) cake
    write( temporary_unit, '("teapot = ''", A, "''")' ) key_from_teapot( teapot )
    write( temporary_unit, '("cheese = ", L2)' ) cheese
    write( temporary_unit, '("fish = ", E14.7)' ) fish
    write( temporary_unit, '("tail = ", I0)' ) tail
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_everything_namelist( temporary_unit, local_rank )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop &
        'feign_everything_config: '// &
        'Unable to close temporary file'

  end subroutine feign_everything_config

end module everything_mod
        '''.strip()

        everything = namelist.NamelistDescription('everything')
        everything.add_string('cake', configure_string_length='filename')
        everything.add_enumeration('teapot', enumerators=['brown', 'steel'])
        everything.add_value('cheese', 'logical')
        everything.add_value('fish', 'real')
        everything.add_usage('wibble', 'constants_mod')
        everything.add_computed('yarn', 'real', 'default',
                                calculation=['fish * wibble / 180.0_r_def'])
        everything.add_value('tail', 'integer', 'native')

        output_file = StringIO.StringIO()
        uut = feigner.NamelistFeigner('everything_mod')
        uut.add_namelist([everything])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())

    ###########################################################################
    def test_multi_file(self):
        expected_source = '''
!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
! Handles the feigning of namelists.
!
module multifile_mod

  use constants_mod, only : i_native, l_def, r_def, str_max_filename
  use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR
  use mpi_mod,       only : get_comm_rank

  implicit none

  private
  public :: feign_first_config, &
            feign_second_config

  integer(i_native) :: local_rank = -1
  integer(i_native), parameter :: temporary_unit = 3

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_first_config( cake, &
                                 teapot, &
                                 cheese )

    use first_config_mod, only : read_first_namelist, &
                                 key_from_teapot, &
                                 teapot_from_key

    implicit none

    character(*), intent(in) :: cake
    integer(i_native), intent(in) :: teapot
    logical(l_def), intent(in) :: cheese

    integer(i_native) :: condition

    if (local_rank == -1) then
      local_rank = get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_first_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&first")' )
    write( temporary_unit, '("cake = ''", A, "''")' ) cake
    write( temporary_unit, '("teapot = ''", A, "''")' ) key_from_teapot( teapot )
    write( temporary_unit, '("cheese = ", L2)' ) cheese
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_first_namelist( temporary_unit, local_rank )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop &
        'feign_first_config: '// &
        'Unable to close temporary file'

  end subroutine feign_first_config

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine feign_second_config( fish, &
                                  yarn, &
                                  tail )

    use second_config_mod, only : read_second_namelist, &
                                  key_from_yarn, &
                                  yarn_from_key

    implicit none

    real(r_def), intent(in) :: fish
    integer(i_native), intent(in) :: yarn
    integer(i_native), intent(in) :: tail

    integer(i_native) :: condition

    if (local_rank == -1) then
      local_rank = get_comm_rank()
    end if

    open( temporary_unit, status='scratch', action='readwrite', &
          iostat=condition )
    if (condition /= 0) then
      write( 6, '("feign_second_config: ", I0)' ) condition
      stop
    end if

    write( temporary_unit, '("&second")' )
    write( temporary_unit, '("fish = ", E14.7)' ) fish
    write( temporary_unit, '("yarn = ''", A, "''")' ) key_from_yarn( yarn )
    write( temporary_unit, '("tail = ", I0)' ) tail
    write( temporary_unit, '("/")' )

    rewind(temporary_unit)
    call read_second_namelist( temporary_unit, local_rank )

    close(temporary_unit, iostat=condition )
    if (condition /= 0) stop &
        'feign_second_config: '// &
        'Unable to close temporary file'

  end subroutine feign_second_config

end module multifile_mod
        '''.strip()

        firstfile = namelist.NamelistDescription('first')
        firstfile.add_string('cake', configure_string_length='filename')
        firstfile.add_enumeration('teapot', enumerators=['brown', 'steel'])
        firstfile.add_value('cheese', 'logical')

        secondfile = namelist.NamelistDescription('second')
        secondfile.add_value('fish', 'real')
        secondfile.add_enumeration('yarn', enumerators=['fuzzy', 'colourful'])
        secondfile.add_value('tail', 'integer', 'native')

        output_file = StringIO.StringIO()
        uut = feigner.NamelistFeigner('multifile_mod')
        uut.add_namelist([firstfile, secondfile])
        uut.write_module(output_file)

        self.assertMultiLineEqual(expected_source + '\n',
                                  output_file.getvalue())
