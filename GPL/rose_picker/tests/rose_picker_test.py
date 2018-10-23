#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# (C) British Crown Copyright 2018 Met Office.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rose. If not, see <http://www.gnu.org/licenses/>.
##############################################################################

import unittest
import os
import collections
import pickle
import json
import tempfile
from subprocess import check_output, Popen, CalledProcessError, STDOUT


PICKER_EXE = './rose_picker'


###############################################################################
class RosePickerTest(unittest.TestCase):

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
    def test_no_namelist_for_member(self):

        input_file = tempfile.NamedTemporaryFile()
        input_file.write('''
[namelist:kevin=orphan]
type=integer
''')
        input_file.seek(0)
        picker_command = "{} {}".format(PICKER_EXE, input_file.name)

        with self.assertRaises(CalledProcessError) as context:
            check_output(picker_command, shell=True, stderr=STDOUT)

        input_file.close()

        self.assertIn(
            'namelist:kevin has no section in metadata configuration file',
            context.exception.output)

    ##########################################################################
    def test_good_picker(self):

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
        input_file.close()

        # json file
        self.nml_config_file = \
            '{}.json'.format(os.path.basename(input_file.name))

        config_file = open(self.nml_config_file, 'rb')
        result = json.load(config_file)
        config_file.close()

        good_result = collections.OrderedDict(
            {'aerial': {'dino':   {'length': ':',
                                   'type':   'integer',
                                   'bounds': 'namelist:sugar=TABLET'},
                        'wilma':  {'length': ':',
                                   'type':   'real',
                                   'bounds': 'source:constants_mod=FUDGE'},
                        'betty':  {'length': ':',
                                   'type':   'logical',
                                   'bounds': 'fred'},
                        'bambam': {'length': ':',
                                   'type':   'integer'},
                        'fred':   {'type':   'real'}}})

        self.assertEqual(good_result, result)
