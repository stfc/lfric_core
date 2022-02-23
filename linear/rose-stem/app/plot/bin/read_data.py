#!/usr/bin/env python
# -*- coding: utf-8 -*-
##############################################################################
# (c) Crown copyright 2022 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
'''
Python script for reading LFRic diagnostic output data in UGRID NetCDF format
'''
import iris


def load_cube_by_varname(filename, var):
    '''
    Read the file and extract field of interest
    '''
    variable_constraint = iris.Constraint(cube_func=(
        lambda c: c.var_name == var))
    return iris.load_cube(filename, constraint=variable_constraint)


def read_ugrid_data(filestem, field):
    '''
    Read from UGRID NetCDF file
    # filestem - path to file
    # field - field name
    '''
    cube = load_cube_by_varname(filestem, field)

    return cube
