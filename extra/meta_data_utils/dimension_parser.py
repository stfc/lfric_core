##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""This module contains all functionality related to vertical dimension meta
data in LFRic. The two main features are the ability to parse the default and
hard coded values of vertical dimensions and to expand vertical dimension
declarations in LFRic meta data
e.g. model_height_dimension(top=TOP_WET_LEVEL)"""
import logging
import re
from typing import Dict, List

from fparser.common.readfortran import FortranFileReader
from fparser.two.Fortran2003 import Assignment_Stmt, Else_Stmt, \
    Execution_Part, Function_Stmt, Function_Subprogram, \
    If_Construct, If_Then_Stmt, Name, Structure_Constructor_2
from fparser.two.parser import ParserFactory
from fparser.two.utils import walk

LOGGER = logging.getLogger(__name__)

DIMENSION_TYPE_REGEX = re.compile(r"(?P<type>[a-zA-Z_]+)[\s]*\([^)]*\)")
TOP_ARG_REGEX = re.compile(r"top[\s=]*(?P<top_arg>[A-Za-z_]+)")
BOTTOM_ARG_REGEX = re.compile(r"bottom[\s=]*(?P<bottom_arg>[A-Za-z_]+)")
LEVEL_DEF_REGEX = re.compile(r"(?P<level>[\d.]+)")
ENUM_REGEX = re.compile(r"(?P<levels>[a-zA-Z_]+)")

F2003_PARSER = ParserFactory().create(std="f2003")


def get_default_values(if_construct, path, levels):
    """Takes an fparser if_construct object as an argument. This if statement
    takes the form of
    if(present(DUMMY_ARG)) then
      a_level = DUMMY_ARG
    else
      a_level = DEFAULT_VALUE
    endif

    This function returns the DUMMY_ARG and DEFAULT_VALUE as a tuple. If values
    not found, the error is logged and a tuple containing none's is returned
    """
    dummy_arg = None
    default_value = None

    for i in range(len(if_construct.content)):

        # Get Dummy argument
        if isinstance(if_construct.content[i], If_Then_Stmt):
            line = if_construct.content[i + 1].items
            dummy_arg = line[2].tostr()

        # Get Default value
        if isinstance(if_construct.content[i], Else_Stmt):
            line = if_construct.content[i + 1].items
            default_value = line[2].tostr()

    if not dummy_arg or not default_value:
        LOGGER.error("File at %s is invalid. Problem with default "
                     "values", path)
    if default_value not in levels:
        LOGGER.error("File at %s is invalid. Default level does not"
                     " exist", path)

    return dummy_arg, default_value


def get_hard_coded_values(assignment_statement):
    """Takes an Fparser assignment_statement object as an argument. It returns
    the object construction arguments as key value pairs.
    The assignment statement takes the form of

    some_var = some_object_constructor(&
                    key = value,                              &
                    key = value)

    :param assignment_statement: An fparser assignment_statement
    :return hard_coded_values: A dictionary containing the keys and values used
    in the assignment statement"""
    hard_coded_values = {}
    unwanted_values = ['upper_level', 'lower_level', 'level_definition']

    for value in walk(assignment_statement, types=Structure_Constructor_2):
        if value.children[0].string in unwanted_values:
            continue

        stripped_value = value.children[1].string.replace("'", "")
        hard_coded_values[value.children[0].string] = stripped_value

    return hard_coded_values


def parse_vertical_dimension_info(path: str, levels: List) -> Dict:
    """Parses though a "vertical_dimension_mod.f90" fortran file. Returns a
    dictionary containing information about functions in the file.
    :param path: Path to the LFRic meta data utility functions
    :param levels: A list of model levels
    :return functions: A dictionary containing vertical dimension definitions
"""

    functions = {}
    reader = FortranFileReader(path + "vertical_dimensions_mod.f90")
    parse_tree = F2003_PARSER(reader)

    # Loop over every function in the file
    for function in walk(parse_tree.content, Function_Subprogram):

        function_dict = {}
        function_name = ""
        default_values = {}

        for function_part in function.children:

            # Get function name
            if isinstance(function_part, Function_Stmt):
                for statement in function_part.items:
                    if isinstance(statement, Name):
                        # Ignore the custom dimensions, no default values
                        if "custom" in statement.tostr():
                            break
                        function_name = statement.tostr()

            elif isinstance(function_part, Execution_Part):

                for statement in function_part.content:

                    # Get hard coded values
                    if isinstance(statement, Assignment_Stmt):
                        function_dict["hard_coded_values"] = \
                            get_hard_coded_values(statement)

                    # Get default values and it's corresponding dummy arg
                    elif isinstance(statement, If_Construct):
                        dummy_arg, default_value = \
                            get_default_values(statement, path, levels)

                        default_values[dummy_arg] = default_value

        function_dict["default_values"] = default_values
        functions[function_name] = function_dict

    return functions


def translate_vertical_dimension(dimension_definition, dimension_declaration):
    """Takes dimension definition as a string and returns a dictionary
    containing that dimension's attributes
    :param dimension_definition: A dictionary containing hard-coded default
    values and for various types of vertical dimension
    :param dimension_declaration: A string that defines the type of vertical
    and any arguments that it is taking
    :return parsed_definition: """

    LOGGER.debug("Parsing a vertical dimension")
    LOGGER.debug("Dimension declaration: %s", dimension_declaration)

    parsed_definition = {}
    dimension_type_match = DIMENSION_TYPE_REGEX.search(dimension_declaration)
    top_arg = TOP_ARG_REGEX.search(dimension_declaration)
    bottom_arg = BOTTOM_ARG_REGEX.search(dimension_declaration)
    fixed_levels = LEVEL_DEF_REGEX.findall(dimension_declaration)

    dimension_type = dimension_type_match.group('type')

    if "model" in dimension_type:
        if top_arg:
            parsed_definition["top_arg"] = top_arg.group("top_arg")
        else:
            parsed_definition["top_arg"] = \
                dimension_definition[dimension_type]["default_values"]["top"]

        if bottom_arg:
            parsed_definition["bottom_arg"] = bottom_arg.group("bottom_arg")
        else:
            parsed_definition["bottom_arg"] = \
                dimension_definition[
                    dimension_type]["default_values"]["bottom"]

    for key, value in dimension_definition[dimension_type][
            "hard_coded_values"].items():

        parsed_definition[key] = value

    if fixed_levels:
        levels = []
        for level in fixed_levels:
            levels.append(float(level))
        parsed_definition["level_definition"] = levels

    LOGGER.debug("Parsed Definition: %s", parsed_definition)
    return parsed_definition
