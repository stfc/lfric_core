#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to strip the crun info.
'''
from jinja2 import contextfilter
import re
import ast

@contextfilter
def getCrunInfo(context, call):
    '''
    Takes a string return dictionary values related to crunning:
        crun: Number of runs to do in the crun.
        ...
    @param [in] context Jinja2 instance to run macro against.
    @return String resulting from setting the environment.
    '''
    if call.find('(') == -1:
        macroName = call
        arguments = ''
    else:
        macroName = call[:call.index('(')]
        arguments = re.split(', *', call[call.index('(')+1:call.rindex(')')])

    normalArguments  = [argument for argument in arguments \
                        if argument.find('=') == -1]
    keywordArguments = [argument for argument in arguments \
                        if argument.find('=') != -1]

    argumentList = []
    for argument in normalArguments:
        if argument[0] == '"':
            argumentList.append( argument[1:-1] )
        else:
            argumentList.append( argument )

    argumentDictionary = {}
    for argument in keywordArguments:
        key, value = re.split(' *= *', argument)
        argumentDictionary[key] = value

    # Return info about the crun arguments
    return_value={}
    if 'crun' in argumentDictionary.keys():
        return_value.update({'crun':int(argumentDictionary['crun'])})

    return return_value
