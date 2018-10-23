#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Turns namelist descriptions into namelist modules.
'''

from __future__ import print_function

from abc import ABCMeta, abstractmethod

import collections
import random
import re
import json

import jinja2 as jinja
import configurator.jinjamacros as jinjamacros


##############################################################################
class NamelistDescriptionException(Exception):
    pass


##############################################################################
class FortranType(object):
    _singletonMap = {}

    @classmethod
    def instance(cls, intrinsic_type, kind, write_format):
        if intrinsic_type not in cls._singletonMap:
            cls._singletonMap[intrinsic_type] = {}

        if kind not in cls._singletonMap[intrinsic_type]:
            cls._singletonMap[intrinsic_type][kind] = {}

        if write_format not in cls._singletonMap[intrinsic_type][kind]:
            cls._singletonMap[intrinsic_type][kind][write_format] \
                = cls(intrinsic_type, kind, write_format)

        return cls._singletonMap[intrinsic_type][kind][write_format]

    def __init__(self, intrinsic_type, kind, write_format):
        self.intrinsic_type = intrinsic_type
        self.kind = kind
        self.write_format = write_format

    def declaration(self):
        return self.intrinsic_type + '(' + self.kind + ')'

    def label(self):
        return self.intrinsic_type + '_' + self.kind

    def __lt__(self, other):
        return self.declaration() < other.declaration()

    def __eq__(self, other):
        self.declaration() == other.declaration()


##############################################################################
class _Property(object):
    __metaclass__ = ABCMeta

    def __init__(self, name, fortran_type):
        self.name = name
        self.fortran_type = fortran_type

    def required_kinds(self):
        '''Returns a list of the fortran "kind" statements used by
           this fortran namelist module'''
        return [self.fortran_type.kind]

    @abstractmethod
    def get_configure_type(self):
        pass


##############################################################################
class _String(_Property):
    _fortranStringMap = {'default':  'str_def',
                         'filename': 'str_max_filename'}

    def __init__(self, name, length=None):

        if not length:
            length = 'default'

        super(_String, self).__init__(
            name, FortranType.instance('character',
                                       self._fortranStringMap[length],
                                       'A'))
        self._missing_data_indicator = 'cmdi'

    def get_configure_type(self):
        return 'string'


##############################################################################
class _Enumeration(_Property):

    def __init__(self, name, keyDictionary):
        super(_Enumeration, self).__init__(name,
                                           FortranType.instance('integer',
                                                                'i_native',
                                                                'I0'))

        self.mapping = keyDictionary
        self.inverse_mapping = {value:
                                key for key, value in self.mapping.iteritems()}
        self.first_key = self.inverse_mapping[min(self.inverse_mapping.keys())]
        self._missing_data_indicator = 'emdi'

    def required_kinds(self):
        return [self.fortran_type.kind, 'str_def']

    def get_configure_type(self):
        return 'enumeration'


##############################################################################
class _Scalar(_Property):
    _fortranKindMap = {'logical': {'default': 'l_def',
                                   'native':  'l_native'},
                       'integer': {'default': 'i_def',
                                   'native':  'i_native',
                                   'short':   'i_short',
                                   'long':    'i_long'},
                       'real':    {'default': 'r_def',
                                   'native':  'r_native',
                                   'single':  'r_single',
                                   'double':  'r_double'}}

    _fortranFormatMap = {'logical': 'L2',
                         'integer': 'I0',
                         'real':    'E14.7'}

    _fortranMissingDataIndicator = {'logical': '.false.',
                                    'integer': 'imdi',
                                    'real':    'rmdi'}

    def __init__(self, name, configure_type, configure_kind=None):

        if not configure_kind:
            configure_kind = 'default'

        super(_Scalar, self).__init__(
            name, FortranType.instance(
                configure_type,
                self._fortranKindMap[configure_type][configure_kind],
                self._fortranFormatMap[configure_type]))
        self._missing_data_indicator = \
            self._fortranMissingDataIndicator[configure_type]

    def get_configure_type(self):
        return 'scalar'

    def get_missing_data_indicator(self):
        return self._missing_data_indicator


##############################################################################
class _Computed(_Scalar):

    def __init__(self, name, configure_type, configure_kind, computation):

        super(_Computed, self).__init__(name, configure_type, configure_kind)
        self.computation = computation[0]

    def get_configure_type(self):
        return 'computed'


##############################################################################
class _Array(_Property):
    def __init__(self, name, contentProperty, bounds):
        super(_Array, self).__init__(name, contentProperty.fortran_type)
        self.content = contentProperty

        if not len(bounds) == 1:
            message = 'Only 1D arrays allowed in configuration: {}'
            raise NamelistDescriptionException(message.format(bounds))

        if ':' in bounds[0] and bounds[0].strip() != ':':
            lower, upper = bounds[0].split(':')

            if (lower.strip() not in ['1', '']):
                message = 'Only lower bound of 1 '\
                          'is allowed in configuration: {}'
                raise NamelistDescriptionException(message.format(bounds[0]))

            self.bounds = upper
        else:
            self.bounds = bounds[0]

    def get_configure_type(self):
        return 'array'

    def is_immdeiate_size(self):
        if self.bounds.isdigit():
            return True

        return False

    def is_deferred_size(self):
        if not self.bounds[0].isdigit() and self.bounds[0] != ':':
            return True

        return False

    def is_arbitrary_size(self):
        if self.bounds[0] == ':':
            return True

        return False


##############################################################################
class NamelistDescription(object):

    def __init__(self, listname):
        self._listname = listname

        self._engine = jinja.Environment(
            loader=jinja.PackageLoader('configurator', 'templates'),
            extensions=['jinja2.ext.do'])
        self._engine.filters['decorate'] = jinjamacros.decorate_macro

        self._parameters = collections.OrderedDict()
        self._module_usage = collections.defaultdict(set)
        self._module_usage['constants_mod'] = set(['cmdi', 'emdi',
                                                   'imdi', 'rmdi'])
        self._enum_pool = list(range(1, 1000))

    def get_namelist_name(self):
        '''Returns the fortran namelist name as string'''
        return self._listname

    def get_module_name(self):
        '''Returns the fortran namelist module name as string'''
        return self._listname + '_config_mod'

    def add_enumeration(self, name, enumerators):
        '''Adds an enumeration variable to the namelist description'''
        if not isinstance(enumerators, list):
            message = 'Expected list of enumerators'
            raise NamelistDescriptionException(message)

        key_dict = collections.OrderedDict()
        for key in enumerators:
            key_dict[key] = random.choice(self._enum_pool)
            self._enum_pool.remove(key_dict[key])

        self._parameters[name] = _Enumeration(name, key_dict)

    def add_usage(self, name, module=None):
        '''Add add variable from a different fortran module that
           this fortran namelist module should access via the fortran
           `use` statement'''
        self._module_usage[module].add(name)

    def add_string(self, name, configure_string_length=None, bounds=None):
        '''Add a scalar/array string variable to the namelist description'''

        new_parameter = _String(name, configure_string_length)

        if bounds:
            bounds[0] = self._dereference_expression(bounds[0])
            self._parameters[name] = _Array(name, new_parameter, bounds)
        else:
            self._parameters[name] = new_parameter

    def add_value(self, name, configure_type, configure_kind=None,
                  bounds=None):
        '''Add a scalar/array variable of type [logical,integer,real]
           to the namelist description'''

        new_parameter = _Scalar(name, configure_type, configure_kind)

        if bounds:
            bounds[0] = self._dereference_expression(bounds[0])
            self._parameters[name] = _Array(name, new_parameter, bounds)
        else:
            self._parameters[name] = new_parameter

    def add_computed(self, name, configure_type, configure_kind=None,
                     calculation=None):
        '''Add a variable the namelist module which is derived from an
           expression provided as a string variable'''

        calculation[0] = self._dereference_expression(calculation[0])
        self._parameters[name] = _Computed(name, configure_type,
                                           configure_kind, calculation)

    def get_parameters(self):
        return self._parameters.values()

    def write_module(self, file_object):

        if not self._parameters:
            message = ('Cannot write a module to load an empty namelist ('
                       + self._listname + ')')
            raise NamelistDescriptionException(message)

        all_kinds = set(['i_native'])
        lone_kind_index = {}
        lone_kind_tally = collections.defaultdict(int)
        namelist = []

        for name, parameter in self._parameters.items():

            all_kinds.update(parameter.required_kinds())

            if not isinstance(parameter, _Computed) and \
               not isinstance(parameter, _Array):

                lone_kind_tally[parameter.fortran_type] += 1
                lone_kind_index[name] = lone_kind_tally[parameter.fortran_type]

            if not isinstance(parameter, _Computed):
                namelist.append(parameter.name)

        inserts = {'all_kinds':     all_kinds,
                   'arrays':        [parameter.name
                                     for parameter in self._parameters.values()
                                     if isinstance(parameter, _Array)],
                   'allocatables':  [parameter.name
                                     for parameter in self._parameters.values()
                                     if (isinstance(parameter, _Array) and
                                         not parameter.is_immdeiate_size())],
                   'enumerations':  [parameter.name
                                     for parameter in self._parameters.values()
                                     if isinstance(parameter, _Enumeration)],
                   'listname':      self._listname,
                   'lonekindindex': lone_kind_index,
                   'lonekindtally': lone_kind_tally,
                   'namelist':      namelist,
                   'parameters':    self._parameters,
                   'use_from':      self._module_usage}

        template = self._engine.get_template('namelist.f90.jinja')
        print(template.render(inserts), file=file_object)

    def _dereference_expression(self, string):

        str_dict = {'namelist': {'regexString':   r'namelist:(\w*)=(\w*)',
                                 'removalString': r'namelist:\w*=',
                                 'moduleSuffix':  '_config_mod'},
                    'source':   {'regexString':   r'source:(\w*)=(\w*)',
                                 'removalString': r'source:\w*=',
                                 'moduleSuffix':  ''}}
        result = string

        for value in str_dict.values():

            use_variables = re.findall(value['regexString'], result)
            if use_variables is not None:
                n_vars = len(use_variables)
                for i_var in range(0, n_vars):
                    if use_variables[i_var][0] != self._listname:
                        module_name = '{}{}'.format(
                            use_variables[i_var][0],
                            value['moduleSuffix'])
                        var_name = use_variables[i_var][1]

                        self.add_usage(var_name, module=module_name)

            result = re.sub(value['removalString'], '', result)

        return result

    def add_member(self, member_name, meta_dict):

        meta_keys = meta_dict.keys()
        string_length = None
        xkind = None
        xtype = None
        xbounds = None

        if 'string_length' in meta_keys:
            string_length = meta_dict['string_length']

        if 'kind' in meta_keys:
            xkind = meta_dict['kind']

        if 'type' in meta_keys:
            xtype = meta_dict['type']
            if isinstance(xtype, str):
                xtype = xtype.replace('character', 'string')

        elif ('enumeration' not in meta_keys or
              meta_dict['enumeration'] == 'false'):
            message = ('namelist:' + self._listname + '='
                       + member_name
                       + ': Non-enumeration metadata requires '
                       + 'a type definition')
            raise NamelistDescriptionException(message)

        # Determining array bounds if any.
        if 'length' in meta_keys:

            xlength = meta_dict['length']

            if xlength == ':':

                if 'bounds' in meta_keys:
                    xbounds = meta_dict['bounds']
                    xbounds = [xbounds]
                else:
                    xbounds = [':']

            elif isinstance(int(xlength), int):
                xbounds = [xlength]

        # Generating Enumerators from metadata
        # These are not dependant on xtype being specified
        if ('enumeration' in meta_keys and
                meta_dict['enumeration'] == 'true'):

            enumeration_keys = meta_dict['values']
            if all(isinstance(item, str) for item in enumeration_keys):
                enumeration_keys = enumeration_keys.replace('\n', '')
                enumeration_keys = enumeration_keys.replace(' ', '')
                enumeration_keys = enumeration_keys.replace("'", '')
                enumeration_keys = enumeration_keys.split(',')

                enumeration_keys = [re.sub(r'namelist:', '', member)
                                    for member in enumeration_keys]

                self.add_enumeration(
                    member_name, enumerators=enumeration_keys)

        # Check to see if member is a derived variable
        elif 'expression' in meta_keys:
            expression_string = meta_dict['expression']
            self.add_computed(
                member_name, xtype, configure_kind=xkind,
                calculation=[expression_string])

        elif xtype == 'string':
            self.add_string(
                member_name,
                configure_string_length=string_length,
                bounds=xbounds)
        else:
            self.add_value(
                member_name, xtype, configure_kind=xkind,
                bounds=xbounds)


###############################################################################
class NamelistConfigDescription(object):

    def __init__(self):
        pass

    def decode_unicode(self, namelist_config):

        if isinstance(namelist_config, dict):
            new_dict = {}
            for key, value in namelist_config.iteritems():
                if isinstance(key, unicode):
                    key = key.encode('ascii')
                if isinstance(value, dict):
                    value = self.decode_unicode(value)
                else:
                    value = value.encode('ascii')

                new_dict[key] = value

            result = new_dict
        else:
            message = 'Supplied namelist configuration ' + \
                      'object is not a dictionary'
            raise NamelistDescriptionException(message)

        return result

    def process_config(self, nml_config_file):

        # Process json
        with open(nml_config_file) as config_file:
            namelist_config = json.load(config_file)

        namelist_config = self.decode_unicode(namelist_config)
        result = []

        for listname in namelist_config.keys():
            description = NamelistDescription(listname)
            list_dict = namelist_config[listname]

            for member in sorted(list_dict.keys()):

                meta_dict = list_dict[member]
                description.add_member(member, meta_dict)

            result.append(description)

        return result
