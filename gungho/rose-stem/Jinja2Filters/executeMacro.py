#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################
'''
Implements a Jinja2 filter to run a macro specified by a string.
'''
from jinja2 import contextfilter
import re

def checkBoolean(string):
    ''' If string is a boolean or None, i.e. "False", "True" or "None", then
    return the appropriate boolean value or None rather than the string '''
    if string=="True": string=True
    if string=="False": string=False
    if string=="None": string=None
    return string

@contextfilter
def executeMacro(context, call):
    '''
    Takes a string and executes it as though it were a Jinja2 macro call

    The call string has the syntax <macro name>([<argument>]...).

    Arguments can be either position or keyword based.

    @param [inout] context Jinja2 instance to run macro against.
    @param [in]    call    Invokation string.
    @return String resulting from calling the macro.
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
            argumentList.append( checkBoolean(argument[1:-2]) )
        else:
            argumentList.append( checkBoolean(argument) )

    argumentDictionary = {}
    for argument in keywordArguments:
        key, value = re.split(' *= *', argument)
        argumentDictionary[key] = checkBoolean(value)

    return context.vars[macroName]( *argumentList, **argumentDictionary )
