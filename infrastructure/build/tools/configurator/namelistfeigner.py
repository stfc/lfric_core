#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################

from __future__ import print_function

from __future__ import absolute_import
import collections
import jinja2 as jinja

import configurator.jinjamacros as jinjamacros


##############################################################################
class NamelistFeigner(object):
    def __init__(self, moduleName):
        self._module_name = moduleName

        self._engine = jinja.Environment(
            loader=jinja.PackageLoader('configurator', 'templates'),
            extensions=['jinja2.ext.do'])
        self._engine.filters['decorate'] = jinjamacros.decorate_macro

        self._namelists = collections.OrderedDict()

    def add_namelist(self, namelists):
        for item in namelists:
            self._namelists[item.get_namelist_name()] = item

    def write_module(self, module_file):
        enumerations = collections.defaultdict(list)
        kinds = set(['i_native'])
        namelists = []
        parameters = {}
        character_arrays = None
        non_character_arrays = None

        for namelist in self._namelists.values():
            namelists.append(namelist.get_namelist_name())
            parameters[namelist.get_namelist_name()] = []
            for param in namelist.get_parameters():
                if param.get_configure_type() == 'enumeration':
                    enumerations[namelist.get_namelist_name()].append(
                        param.name)
                if param.get_configure_type() != 'computed':
                    parameters[namelist.get_namelist_name()].append(param)
                    kinds.add(param.fortran_type.kind)
                if param.get_configure_type() == 'array':
                    if param.fortran_type.intrinsic_type == 'character':
                        character_arrays = True
                    else:
                        non_character_arrays = True

        inserts = {'enumerations': enumerations,
                   'kinds':        kinds,
                   'modulename':   self._module_name,
                   'namelists':    namelists,
                   'parameters':        parameters,
                   'string_arrays':     character_arrays,
                   'non_string_arrays': non_character_arrays}

        template = self._engine.get_template('feign_config.f90.jinja')
        print(template.render(inserts), file=module_file)
