#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
# Process previously analysed dependency database. For fun and profit!

from __future__ import absolute_import, print_function

import logging
import os.path

from utilities.path import replaceExtension

###############################################################################
# Process dependency database.
#
class FortranProcessor():
    ###########################################################################
    # Constructor.
    #
    # Arguments:
    #   database  - FortranDatabase object holding details.
    #   objectdir - The directory which holds .o files.
    #   moduledir - The directory which holds .mod files.
    #
    def __init__(self, database, objectDirectory, moduleDirectory):
        self._database        = database
        self._objectDirectory = objectDirectory
        self._moduleDirectory = moduleDirectory

    ###########################################################################
    # Examine the program unit dependecies and work out the file dependencies.
    #
    # Arguments:
    #   fileStore - FileDependencies object to accept computed dependencies.
    #
    def determineCompileDependencies( self, fileStore ):
        for unit, unitFilename, prerequisite, prerequisiteFilename \
                in self._database.getCompileDependencies():
            message = '{0} depends on {1}'.format(unit, prerequisite)
            logging.getLogger(__name__).info(message)

            objectFilename = replaceExtension( unitFilename, 'o' )
            objectPathname = os.path.join( self._objectDirectory, \
                                           objectFilename )

            moduleFilename = replaceExtension( prerequisiteFilename, 'mod' )
            modulePathname = os.path.join( self._moduleDirectory, \
                                           moduleFilename )

            fileStore.addFileDependency( objectPathname, modulePathname )

    ###########################################################################
    # Determine all program units needed to build each program.
    #
    # TODO: Once we have a more recent version of SQLite we could consider
    # doing this at the database level.
    #
    def determineLinkDependencies( self ):
        for program in self._database.getPrograms():
            logging.getLogger(__name__).info('Program {0}'.format(program))

            prerequisites = set()
            self._descend( program, prerequisites )
            unit, unit_file, prereq, prereq_file = next(self._database.getLinkDependencies( program ))
            program_object_file = os.path.join( self._objectDirectory, \
                                          replaceExtension( unit_file, 'o' ) )
            prerequisites.add( program_object_file )
            yield os.path.join( self._objectDirectory, program ), \
                  list(prerequisites)

    ##########################################################################
    def _descend( self, programUnit, prerequisites ):
        logger = logging.getLogger(__name__)
        for unit, unit_file, prereq, prereq_file \
                         in self._database.getLinkDependencies( programUnit ):
            logger.info('  Requires {0}'.format(unit))

            prereq_object_file = os.path.join( self._objectDirectory, \
                                        replaceExtension( prereq_file, 'o' ) )

            if prereq_object_file in prerequisites:
                logger.info('    Seen already, stopping descent')
            else:
                prerequisites.add( prereq_object_file )
                self._descend( prereq, prerequisites )
