#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from __future__ import print_function

from __future__ import absolute_import
import os.path
import shutil
import tempfile
import unittest

import dependerator.analyser
import dependerator.database
import six

##############################################################################
class NamelistDescriptionAnalyserTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        dbFilename = os.path.join( self._scratchDirectory, 'test.db' )
        self._database     = dependerator.database.SQLiteDatabase( dbFilename )
        self._dependencies = dependerator.database.FileDependencies( self._database )

    ##########################################################################
    def tearDown( self ):
        del self._dependencies
        del self._database
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    # Presents a namelist description to the analyser and ensures the correct
    # dependencies are returned.
    #
    def testAnalysis( self ):
        testFilename = os.path.join( self._scratchDirectory, 'test.nld' )
        with open( testFilename, 'wt' ) as nldFile:
            print( '''
namelist foo

  bar  : string( filename ) !gumph
  baz  : enumeration[ thing1, thing2 ]
  qux  : real
  fred : real[ 'qux * 2' ]

end namelist foo
                   '''.strip(), file=nldFile )

        uut = dependerator.analyser \
              .NamelistDescriptionAnalyser(self._dependencies)
        uut.analyse( testFilename )
        dependencies = list(self._dependencies.getDependencies())
        self.assertEqual( [(u'foo_configuration_mod.f90',
                           [testFilename])], dependencies )

##############################################################################
class FortranAnalyserTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self.maxDiff = None

        self._scratchDirectory = tempfile.mkdtemp()
        dbFilename = os.path.join( self._scratchDirectory, 'test.db' )
        self._database     = dependerator.database.SQLiteDatabase( dbFilename )
        self._dependencies = dependerator.database \
                             .FortranDependencies( self._database )

    ##########################################################################
    def tearDown( self ):
        del self._dependencies
        del self._database
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    # Ensure continuation lines are handled correctly.
    #
    def testContinuationLines( self ):
      self._dependencies.addProgram( u'stock', u'oxo.f90' )
      self._dependencies.addModuleCompileDependency( u'stock', u'widgit' )
      self._dependencies.addModuleLinkDependency( u'stock', u'widgit' )
      self._dependencies.addModuleCompileDependency( u'stock', u'thingy' )
      self._dependencies.addModuleLinkDependency( u'stock', u'thingy' )
      self._dependencies.addModule( u'beef', u'cow.f90' )
      self._dependencies.addModule( u'pork', u'pig.f90' )
      testFilename = six.text_type( os.path.join( self._scratchDirectory,
                                            'cont.f90' ) )
      with open( testFilename, 'wt' ) as fortranFile:
        print( '''
subroutine thingy(cheese, meat, &
                  teapot, fishslice, &
                 )

end subroutine thingy

subroutine widgit( grunk, &
! This comment should be ignored
                   grank )
  use beef
  implicit none
end subroutine widgit

subroutine old_school( bibble, &
                       & bobble )
  use pork
  implicit none
end subroutine old_school
               '''.strip(), file=fortranFile )
      uut = dependerator.analyser.FortranAnalyser([], self._dependencies)
      uut.analyse( testFilename )

      dependencies = list(self._dependencies.getCompileDependencies())
      self.assertEqual( [(u'old_school', testFilename, u'pork', u'pig.f90'),
                         (u'stock', u'oxo.f90', u'thingy', testFilename),
                         (u'stock', u'oxo.f90', u'widgit', testFilename),
                         (u'widgit', testFilename, u'beef', u'cow.f90')],
                        sorted(dependencies) )

      dependencies = list(self._dependencies.getLinkDependencies( 'widgit' ))
      self.assertEqual( [('widgit', testFilename, u'beef', u'cow.f90')],
                        dependencies )

    ##########################################################################
    # Procedure as program unit.
    #
    def testProcedure( self ):
      testFilename = six.text_type( os.path.join( self._scratchDirectory,
                              'test.f90' ) )
      with open( testFilename, 'wt' ) as fortranFile:
          print( '''
subroutine empty_sub()
end subroutine empty_sub

subroutine one_sub( argument )
end subroutine one_sub

real function cosd(degree)
end function cosd
                  '''.strip(), file=fortranFile )

      uut = dependerator.analyser.FortranAnalyser([], self._dependencies)
      uut.analyse( testFilename )

      self.assertEqual( [(u'cosd',      testFilename),
                         (u'empty_sub', testFilename),
                         (u'one_sub',   testFilename)],
                        sorted(self._dependencies.get_program_units()) )

    ##########################################################################
    # Includes disparate case to ensure case insensitivity.
    #
    def testAnalyseProgram( self ):
        self._dependencies.addModule( 'constants_mod', 'constants_mod.f90' )
        self._dependencies.addModule( 'trumpton_mod', 'trumpton_mod.f90' )

        testFilename = os.path.join( self._scratchDirectory, 'test.f90' )
        with open( testFilename, 'wt' ) as fortranFile:
            print( '''
program fOo

  use constAnts_mod, only : i_def
  use trumpton_Mod, only : hew, pew, barney, mcgrey, cuthbirt, dibble, grub

  implicit none

end program fOo
                   '''.strip(), file=fortranFile )

        uut = dependerator.analyser.FortranAnalyser([], self._dependencies)
        uut.analyse( testFilename )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [u'foo'], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'foo', six.text_type(testFilename),
                            u'constants_mod', u'constants_mod.f90'),
                           (u'foo', six.text_type(testFilename),
                            u'trumpton_mod', u'trumpton_mod.f90')],
                          sorted(dependencies) )

        dependencies = list(self._dependencies.getLinkDependencies( 'foo' ))
        self.assertEqual( [(u'foo', six.text_type(testFilename),
                            u'constants_mod', u'constants_mod.f90'),
                           (u'foo', six.text_type(testFilename),
                            u'trumpton_mod', u'trumpton_mod.f90')],
                          sorted(dependencies) )

    ##########################################################################
    # Includes disparate case to ensure case insensitivity.
    #
    def testAnalyseModule( self ):
        uut = dependerator.analyser.FortranAnalyser([], self._dependencies)

        testFilename = os.path.join( self._scratchDirectory, 'test.f90' )
        with open( testFilename, 'wt' ) as fortranFile:
            print( '''
module foO

  ! Ignore this
  use consTants_mod, only : i_def
  use trumPton_mod, only : hew, pew, barney, mcgrey, cuthbirt, dibble, grub

  implicit none

  private

contains

end module foO

module truMpton_mod

end module truMpton_mod
                   '''.strip(), file=fortranFile )
        uut.analyse( testFilename )

        otherFilename = os.path.join( self._scratchDirectory, 'other.f90' )
        with open( otherFilename, 'wt' ) as otherFile:
            print( '''
module coNstants_mod

end module coNstants_mod
                   '''.strip(), file=otherFile )
        uut.analyse( otherFilename )


        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'foo', testFilename,
                            u'constants_mod', otherFilename),
                           (u'foo', testFilename,
                            u'trumpton_mod', testFilename)],
                          sorted(dependencies) )

        dependencies = list(self._dependencies.getLinkDependencies( 'foo' ))
        self.assertEqual( [(u'foo', testFilename,
                            u'constants_mod', otherFilename),
                           (u'foo', testFilename,
                            u'trumpton_mod', testFilename)],
                          sorted(dependencies) )

   ##########################################################################
    # This test also includes disparate case to ensure case insensitivity is
    # enforced.
    #
    def testAnalyseSubModule( self ):
        uut = dependerator.analyser.FortranAnalyser([], self._dependencies)

        parentFilename = six.text_type(os.path.join( self._scratchDirectory,
                                               'parent.f90' ), "utf-8")
        with open( parentFilename, 'wt' ) as parentFile:
            print( '''
module Parent

  implicit none

  private

  type, public :: test_type
  contains
    procedure foo
    procedure bar
    procedure baz
  end type test_type

  interface
    module subroutine foo( this, cheese )
      class(test_type), intent(inout) :: this
      real,             intent(in)    :: cheese
    end subroutine foo

    module subroutine bar( this, teapot )
      class(test_type), intent(inout) :: this
      character(*),     intent(in)    :: teapot
    end subroutine bar

    type(something) module function baz( this )
      class(test_type), intent(in) :: this
    end function baz
  end interface

end module Parent

submodule (pArent) chIld3

  implicit none

contains

  type(something) module function baz( this )
    class(test_type), intent(in) :: this
  end function baz

end submodule chIld3
                   '''.strip(), file=parentFile )
        uut.analyse( parentFilename )

        child1Filename = six.text_type(os.path.join( self._scratchDirectory,
                                               'child1.f90' ), 'utf-8')
        with open( child1Filename, 'wt' ) as child1File:
            print( '''
submodule (paRent) Child1

  implicit none

  type :: secondary_type
  contains
    procedure baz
  end type secondary_type

  interface
    module subroutine baz( this, wibble )
      class(secondary_type), intent(inout) :: this
      real,                  intent(in)    :: wibble
    end subroutine baz
  end interface

contains

  module subroutine foo( this, cheese )

    implicit none

    class(test_type), intent(inout) :: this
    real,             intent(in)    :: cheese

    type(secondary_type) :: thang

    thang = secondary_type()

    write(6, *) cheese
    call thang%baz( 12.7 )

  end subroutine foo

end submodule Child1
                   '''.strip(), file=child1File )
        uut.analyse( child1Filename )

        child2Filename = six.text_type(os.path.join( self._scratchDirectory,
                                 'child2.f90' ), 'utf-8')
        with open( child2Filename, 'wt' ) as child2File:
            print( '''
submodule (parent) cHild2

  implicit none

contains

  module procedure bar

    implicit none

    write( 6, *) teapot

  end procedure bar

end submodule cHild2
                   '''.strip(), file=child2File )
        uut.analyse( child2Filename )

        child3Filename = six.text_type(os.path.join( self._scratchDirectory,
                                 'child3.f90' ), 'utf-8')
        with open( child3Filename, 'wt' ) as child3File:
            print( '''
submodule (parEnt:chilD1) grandChild

  implicit none

contains

  module subroutine baz( this, wibble )

    implicit none

    class(secondary_type), intent(inout) :: this
    real,                  intent(in)    :: wibble

    write(6, *) wibble

  end subroutine baz

end submodule grandChild
                   '''.strip(), file=child3File )
        uut.analyse( child3Filename )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'child1', child1Filename,
                            u'parent', parentFilename),
                           (u'child2', child2Filename,
                            u'parent', parentFilename),
                           (u'grandchild', child3Filename,
                            u'child1', child1Filename)],
                          sorted(dependencies) )

        dependencies = list(self._dependencies.getLinkDependencies( 'parent' ))
        self.assertEqual( [(u'parent', parentFilename,
                            u'child1', child1Filename),
                           (u'parent', parentFilename,
                            u'child2', child2Filename)],
                          sorted(dependencies) )
        dependencies = list(self._dependencies.getLinkDependencies( 'child1' ))
        self.assertEqual( [(u'child1', child1Filename,
                            u'grandchild', child3Filename)],
                          sorted(dependencies) )

    ##########################################################################
    # Ensure the analyser isn't tripped up by naked global level procedures as
    # program units.
    #
    def testFunctionInModuleName( self ):
        uut = dependerator.analyser.FortranAnalyser([], self._dependencies)

        testFilename = os.path.join( self._scratchDirectory, 'test.f90' )
        with open( testFilename, 'wt' ) as fortranFile:
            print( '''
module function_thing_mod

  use constants_mod, only : i_def

  implicit none

  private

contains

end module function_thing_mod
                   '''.strip(), file=fortranFile )
        uut.analyse( testFilename )

        otherFilename = os.path.join( self._scratchDirectory, 'other.f90' )
        with open( otherFilename, 'wt' ) as otherFile:
            print( '''
module constants_mod

end module constants_mod
                   '''.strip(), file=otherFile )
        uut.analyse( otherFilename )

        dependFilename = os.path.join( self._scratchDirectory, 'dependson.f90' )
        with open( dependFilename, 'wt' ) as dependFile:
            print( '''
subroutine dependson

end dependson
                   '''.strip(), file=dependFile )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'function_thing_mod', testFilename,
                            u'constants_mod', otherFilename)],
                          sorted(dependencies ) )

        dependencies = list(self._dependencies \
                                .getLinkDependencies( 'function_thing_mod' ))
        self.assertEqual( [(u'function_thing_mod', testFilename,
                            u'constants_mod', otherFilename)],
                          sorted(dependencies ) )

    ##########################################################################
    # The analyser has to be able to track dependencies using the deprecated
    # "depends on:" comments of the UM.
    #
    def testDependsOn( self ):
        self._dependencies.add_procedure( u'flibble', u'flibble.f90' )
        testFilename = six.text_type( os.path.join( self._scratchDirectory,
                                              'test.f90' ) )
        with open( testFilename, 'wt' ) as fortranFile:
            print( '''
module function_thing_mod

  use constants_mod, only : i_def

  implicit none

! Add in an interface block - this will test to make sure
! we don't pick up a spurious subroutine call
  interface
     subroutine wooble ()
     end subroutine
  end interface

  private

! Comments before the "depends on" shouldn't upset it.

! depends on: flibble.o

! depends on: wooble

contains

end module function_thing_mod
                   '''.strip(), file=fortranFile )

        otherFilename = six.text_type( os.path.join( self._scratchDirectory,
                                 'other.f90' ) )
        with open( otherFilename, 'wt' ) as otherFile:
            print( '''
module constants_mod
contains
subroutine wooble

end subroutine wooble
end module constants_mod
                   '''.strip(), file=otherFile )

        dependFilename = six.text_type( os.path.join( self._scratchDirectory,
                                  'wooble.f90' ) )
        with open( dependFilename, 'wt' ) as dependFile:
            print( '''
subroutine wooble

end subroutine wooble
                   '''.strip(), file=dependFile )

        uut = dependerator.analyser.FortranAnalyser([], self._dependencies)
        uut.analyse( testFilename )
        uut.analyse( otherFilename )
        uut.analyse( dependFilename )

        programs = list( self._dependencies.getPrograms() )
        self.assertEqual( [], programs )

        dependencies = list(self._dependencies.getCompileDependencies())
        self.assertEqual( [(u'function_thing_mod', testFilename,
                            u'constants_mod', otherFilename),
                           (u'function_thing_mod', testFilename,
                            u'flibble', u'flibble.f90'),
                           (u'function_thing_mod', testFilename,
                            u'wooble', dependFilename)],
                          sorted(dependencies) )

        dependencies = list(self._dependencies \
                                .getLinkDependencies( 'function_thing_mod' ))
        self.assertEqual( [(u'function_thing_mod', testFilename,
                            u'constants_mod', otherFilename),
                           (u'function_thing_mod', testFilename,
                            u'flibble', u'flibble.f90'),
                           (u'function_thing_mod', testFilename,
                            u'wooble', dependFilename)],
                          sorted(dependencies)  )

    ##########################################################################
    # The analyser must ignore abstract interface definitions. These are not
    # program units.
    #
    def testAbstractInterface( self ):
      firstFilename = os.path.join( self._scratchDirectory, 'test.f90' )
      with open( firstFilename, 'wt' ) as fortranFile:
        print( '''
module first_mod

implicit none

private

abstract interface
    subroutine thing_face()
      implicit none
    end subroutine thing_face
end interface

contains

end module first_mod
               '''.strip(), file=fortranFile )

      secondFilename = os.path.join( self._scratchDirectory, 'test2.f90' )
      with open( secondFilename, 'wt' ) as fortranFile:
        print( '''
module second_mod

implicit none

private

abstract interface
    subroutine thing_face()
      implicit none
    end subroutine thing_face
end interface

contains

end module second_mod
               '''.strip(), file=fortranFile )

      uut = dependerator.analyser.FortranAnalyser([], self._dependencies)
      uut.analyse( firstFilename )
      uut.analyse( secondFilename )

    ##########################################################################
    # Ensure "external" works as a dependency marker.
    #
    def testExternal( self ):
      self._dependencies.addModule( 'wibble', 'wibble.f90' )
      self._dependencies.addModule( 'bibble', 'bibble.f90' )
      self._dependencies.addModule( 'ibble', 'ibble.f90' )
      self._dependencies.addModule( 'gribble', 'gribble.f90' )
      testFilename = six.text_type( os.path.join( self._scratchDirectory,
                                            'test.f90' ) )
      with open( testFilename, 'wt' ) as fortranFile:
        print( '''
program boo

  implicit none

  external ibble
  external wibble, bibble, gribble

  call wibble()
  call bibble()

end program boo
               '''.strip(), file=fortranFile )

      uut = dependerator.analyser.FortranAnalyser([], self._dependencies)
      uut.analyse( testFilename )

      programs = list( self._dependencies.getPrograms() )
      self.assertEqual( [u'boo'], programs )

      dependencies = list(self._dependencies.getCompileDependencies())
      self.assertEqual( [(u'boo', testFilename,
                          u'bibble', u'bibble.f90'),
                         (u'boo', testFilename,
                          u'gribble', u'gribble.f90'),
                         (u'boo', testFilename,
                          u'ibble', u'ibble.f90'),
                         (u'boo', testFilename,
                          u'wibble', u'wibble.f90')], sorted(dependencies) )

      dependencies = list(self._dependencies \
                              .getLinkDependencies( 'boo' ))
      self.assertEqual( [(u'boo', testFilename,
                          u'bibble', u'bibble.f90'),
                         (u'boo', testFilename,
                          u'gribble', u'gribble.f90'),
                         (u'boo', testFilename,
                          u'ibble', u'ibble.f90'),
                         (u'boo', testFilename,
                          u'wibble', u'wibble.f90')], sorted(dependencies) )

if __name__ == '__main__':
    unittest.main()
