#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
# For further details please refer to the file LICENCE.original which you
# should have received as part of this distribution.
##############################################################################


###############################################################################
def decorate_macro(subject, prefix=None, postfix=None):
    result = [value for value in subject]

    if prefix:
        result = [prefix+value for value in result]

    if postfix:
        result = [value+postfix for value in result]

    return result
