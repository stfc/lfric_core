###############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
"""Test the dimension parsing..."""
import os

from diag_meta_gen import get_root_dir
from dimension_parser import *
from fortran_reader import read_enum

TEST_DIR = os.path.dirname(os.path.abspath(__file__)) + "/test_data"
ROOT_DIR = get_root_dir() + "/um_physics/source/diagnostics_meta/meta_types/"
ENUM_TEST_FILE = "/enum_test_file"
DIMENSION_TEST_FILE = "/dimension_parser_test_file"
DEFAULT_TEST_FILE = "/get_default_values_test_file"
LEVELS = read_enum(ROOT_DIR + "levels_enum_mod.f90")


def test_read_enum():
    assert read_enum(TEST_DIR + ENUM_TEST_FILE) == \
           ["ONE", "TWO"]


def test_get_default_values():

    reader = FortranFileReader(TEST_DIR + DEFAULT_TEST_FILE)
    parse_tree = F2003_PARSER(reader)

    results = []

    for statement in walk(parse_tree.content, If_Construct):
        results.append(get_default_values(statement, ROOT_DIR, LEVELS))

    assert results == [('test_value_1', 'SOME_TEST_VALUE'),
                       ('Test_value_2', 'SOME_OTHER_TEST_VALUE')]


def test_get_hard_coded_values():
    test_file = "/get_hard_coded_value_test_file"
    reader = FortranFileReader(TEST_DIR + test_file)

    parse_tree = F2003_PARSER(reader)

    results = []

    for statement in walk(parse_tree.content, Assignment_Stmt):
        results.append(get_hard_coded_values(statement))
    assert results == [{
                           'another_test_string': '',
                           'test_name': 'test_string',
                           'test_symbol': 'TEST_SYMBOL'
                           }]


def test_parse_vertical_dimension_info():

    dimension_def = parse_vertical_dimension_info(ROOT_DIR, LEVELS)

    assert "model_depth_dimension" in dimension_def.keys()
    assert "model_height_dimension" in dimension_def.keys()
    assert "fixed_depth_dimension" in dimension_def.keys()
    assert "fixed_height_dimension" in dimension_def.keys()


def test_translate_vertical_dimension():

    dimension_def = parse_vertical_dimension_info(ROOT_DIR, LEVELS)

    model_height = translate_vertical_dimension(dimension_def,
                                                "model_height_dimension("
                                                "bottom="
                                                "BOTTOM_ATMOSPHERIC_LEVEL, "
                                                "top=TOP_ATMOSPHERIC_LEVEL)")
    model_depth = translate_vertical_dimension(dimension_def,
                                               "model_depth_dimension("
                                               "bottom=BOTTOM_SOIL_LEVEL, "
                                               "top=TOP_SOIL_LEVEL)")
    fixed_height = translate_vertical_dimension(dimension_def,
                                                "fixed_height_dimension()")
    fixed_depth = translate_vertical_dimension(dimension_def,
                                               "fixed_depth_dimension()")

    assert model_height == {"standard_name": "height",
                            "units": "m",
                            "top_arg": "TOP_ATMOSPHERIC_LEVEL",
                            "bottom_arg": "BOTTOM_ATMOSPHERIC_LEVEL",
                            "positive": "POSITIVE_UP"}
    assert model_depth == {"standard_name": "depth",
                           "units": "m",
                           "top_arg": "TOP_SOIL_LEVEL",
                           "bottom_arg": "BOTTOM_SOIL_LEVEL",
                           "positive": "POSITIVE_DOWN"}
    assert fixed_height == {"standard_name": "height",
                            "units": "m",
                            "positive": "POSITIVE_UP"}
    assert fixed_depth == {"standard_name": "depth",
                           "units": "m",
                           "positive": "POSITIVE_DOWN"}
