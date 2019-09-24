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

import dependerator.database
import dependerator.process

##############################################################################
class NamelistDescriptionAnalyserTest(unittest.TestCase):
    ##########################################################################
    def setUp( self ):
        self._scratchDirectory = tempfile.mkdtemp()
        dbFilename = os.path.join( self._scratchDirectory, 'fortran.db' )
        database = dependerator.database.SQLiteDatabase( dbFilename )
        self._fortranDatabase = dependerator.database.FortranDependencies( database )
        self._fileDatabase = dependerator.database.FileDependencies( database )

    ##########################################################################
    def tearDown( self ):
        del self._fileDatabase
        del self._fortranDatabase
        shutil.rmtree( self._scratchDirectory )

    ##########################################################################
    def testCompileDependencies( self ):
        self.populateDatabase()

        uut = dependerator.process.FortranProcessor(self._fortranDatabase,
                                                    "objects", "modules")
        uut.determineCompileDependencies( self._fileDatabase )

        self.assertEqual( [(u'objects/bits/bar.o', [u'modules/bits/baz.mod']), \
                           (u'objects/bobs/qux.o', [u'modules/bits/baz.mod']), \
                           (u'objects/foo.o', [u'modules/bits/bar.mod'])],     \
                          list(self._fileDatabase.getDependencies()))

    ##########################################################################
    def testLinkDependencies( self ):
        self.populateDatabase()

        import sys
        uut = dependerator.process.FortranProcessor(self._fortranDatabase,
                                                    "objects", "modules")
        result = list(uut.determineLinkDependencies())

        self.assertEqual( [(u'objects/foo', ['objects/bobs/qux.o', \
                                             'objects/bits/baz.o', \
                                             'objects/bits/bar.o', \
                                             'objects/foo.o'])],   \
                          result )

    ##########################################################################
    def populateDatabase( self ):
        self._fortranDatabase.addProgram( 'foo', 'foo.f90' )
        self._fortranDatabase.addModule( 'bar', 'bits/bar.f90' )
        self._fortranDatabase.addModule( 'baz', 'bits/baz.f90' )
        self._fortranDatabase.addModule( 'qux', 'bobs/qux.f90' )

        self._fortranDatabase.addModuleCompileDependency( 'foo', 'bar' )
        self._fortranDatabase.addModuleCompileDependency( 'bar', 'baz' )
        self._fortranDatabase.addModuleCompileDependency( 'qux', 'baz' )

        self._fortranDatabase.addModuleLinkDependency( 'foo', 'bar' )
        self._fortranDatabase.addModuleLinkDependency( 'bar', 'baz' )
        self._fortranDatabase.addModuleLinkDependency( 'baz', 'qux' )
