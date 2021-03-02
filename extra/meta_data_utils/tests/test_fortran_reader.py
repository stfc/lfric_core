###############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
"""test the fortran reader processes files"""
import os

from diag_meta_gen import get_root_dir
from fortran_reader import *


TEST_DATA_DIR = os.path.dirname(os.path.abspath(__file__)) + "/test_data"
meta_types_folder = "/um_physics/source/diagnostics_meta/meta_types/"
root_dir = get_root_dir()


def test_read_fortran_files_1():
    """does it find the files in the overall project?"""
    test_parser = FortranMetaDataReader(root_dir, meta_types_folder)
    result = test_parser.read_fortran_files()
    assert "example_science_section" in result[0]["sections"].keys()
    # assert result[1] is True # can't test this as it's outside of your
    # control
    assert isinstance(result[0]["sections"]["example_science_section"],
                      Section)


def test_read_fortran_files_2():
    """does our test file from test_data load?"""
    test_parser = FortranMetaDataReader(root_dir, meta_types_folder)
    test_parser.meta_mod_files = [TEST_DATA_DIR +
                                  "/test_section__test_group__.f90"]
    result = test_parser.read_fortran_files()

    assert result[1] is True
    assert "test_section" in result[0]["sections"].keys()
    assert "test_group" in result[0]["sections"]["test_section"].groups
    assert "example_fields__eastward_wind" in \
           result[0]["sections"]["test_section"].groups["test_group"].fields
