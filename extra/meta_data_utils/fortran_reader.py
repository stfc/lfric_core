##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""This module contains functionality for finding and parsing LFRic diagnostic
 meta data files"""
import glob
import logging
import os
import re
from typing import Dict, Tuple

from fparser.common.readfortran import FortranFileReader
from fparser.two.Fortran2003 import Array_Constructor, Assignment_Stmt, \
    Char_Literal_Constant, Enumerator_Def_Stmt, \
    Level_3_Expr, Part_Ref, Section_Subscript_List, Structure_Constructor, \
    Structure_Constructor_2
from fparser.two.parser import ParserFactory
from fparser.two.utils import FparserException, walk

from dimension_parser import parse_vertical_dimension_info, \
    translate_vertical_dimension
from entities import Field, Group, Section
from field_validator import validate_field

F2003_PARSER = ParserFactory().create(std="f2003")


class FortranMetaDataReader:
    """This encapsulates all the parsing functionality. It is give a root
    directory (Head of the repository) upon creation"""

    LOGGER = logging.getLogger("lfric_meta_data_parser")
    vertical_dimension_definition = None
    valid_meta_data = True

    FILE_NAME_REGEX = re.compile(r"(?P<section_name>[a-z_]+?)__"
                                 r"(?P<group_name>[a-z_]+)__")

    def __init__(self, root_directory: str, meta_types_path: str):
        self.__root_dir = root_directory
        self.meta_types_path = meta_types_path
        self.levels = read_enum(self.__root_dir + meta_types_path +
                                "levels_enum_mod.f90")
        self.meta_mod_files = None
        self.find_fortran_files()
        self.get_vertical_dimension_definition()

    def find_fortran_files(self):
        """Recursively looks for fortran files ending with "__meta_mod.f90"
        Initialises the vertical_dimension_definition variable
        Returns a list of file names"""

        self.LOGGER.info("Scanning for fortran meta data files...")
        self.meta_mod_files = glob.glob(
                self.__root_dir +
                '/**/source/diagnostics_meta/**/*__meta_mod.*90',
                recursive=True)

        self.meta_mod_files.sort()

        for file in self.meta_mod_files:
            self.LOGGER.debug("Found meta date file at: " + file)

        self.LOGGER.info("Found %i meta data files", len(self.meta_mod_files))

    def get_vertical_dimension_definition(self):
        """Looks for functions that create vertical dimension objects and reads
        the hard coded values"""
        self.vertical_dimension_definition = \
            parse_vertical_dimension_info(
                    self.__root_dir + self.meta_types_path, self.levels)

    def read_fortran_files(self) -> Tuple[Dict, bool]:
        """Takes a list of file names (meta_mod.f90 files)
        Checks for correctness and returns the relevant fortran lines in a list
        :return Metadata: A dictionary, each key represents a fortran file and
        it's value is a list of strings, each element representing a field"""

        sections_dict: Dict[str, Section] = {}
        valid_files = 0
        # Loop over each found fortran file
        for file_path in self.meta_mod_files:

            try:
                # Load the fortran file

                reader = FortranFileReader(file_path, ignore_comments=True)

                parse_tree = F2003_PARSER(reader)

                file_valid = True
                file_name_parts = self.FILE_NAME_REGEX.search(file_path)

                if not file_name_parts:
                    self.LOGGER.error('Filename in path is not correct' +
                                      os.linesep + file_path)
                    self.valid_meta_data = False
                    break

                section_name = file_name_parts.group("section_name")
                group_name = file_name_parts.group("group_name")
                file_name = file_path[file_path.rfind("/") + 1:]

                if section_name not in sections_dict:
                    sections_dict.update(
                            {section_name: Section(name=section_name)}
                            )

                group = Group(name=group_name, file_name=file_name)

                sections_dict[section_name].add_group(group)

                # For every instance of a meta type object being created
                for definition in walk(parse_tree.content,
                                       types=Assignment_Stmt):

                    field, valid = self.extract_field(definition, file_name)

                    if not valid:
                        self.valid_meta_data = False
                        file_valid = False

                    if validate_field(field):
                        group.add_field(field)
                    else:
                        self.LOGGER.error("%s is invalid. Please check",
                                          field.unique_id)
                        self.valid_meta_data = False
                        file_valid = False

                if file_valid:
                    valid_files += 1

            except FparserException as error:
                self.LOGGER.error(": Fparser Exception %s", error)
                self.valid_meta_data = False

        meta_data = {"sections": sections_dict,
                     "standard_level_markers": self.levels}

        if valid_files == len(self.meta_mod_files):
            self.LOGGER.info("All %i files are valid",
                             len(self.meta_mod_files))
        else:
            self.LOGGER.error("%i of %i files are invalid",
                              len(self.meta_mod_files) - valid_files,
                              len(self.meta_mod_files))

        return meta_data, self.valid_meta_data

    def extract_field(self, definition, file_name) -> Tuple[Field, bool]:
        """Takes an fparser object and extracts the field definition
        information
        :param definition: An fparser object that contains a field definition
        :param file_name: The name of the file that the field is declared in.
        This is needed for Field object creation
        :return Field"""

        valid_field = True
        field = Field(file_name)

        # For every instance argument in object creation
        for parameter in walk(definition, Structure_Constructor_2):
            key = parameter.children[0].string

            # This ignore's args used in vertical dimension creation
            key_blacklist = ["top", "bottom"]
            if key in key_blacklist:
                continue

            # Adds the key / value to the Field object
            if not hasattr(field, key):
                self.logger.error("Unexpected Field Property: %s", key)
                valid_field = False

            try:
                # For multi line statements - override the key value
                if isinstance(parameter.parent, Level_3_Expr):
                    key, value = self.extract_multi_line_statement(parameter)
                    field.add_value(key, value)

                # ENUM's
                elif isinstance(parameter.children[1], Char_Literal_Constant):
                    field.add_value(key, parameter.children[1].string[1:-1])

                # Dimension object creation without args
                # (Structure_Constructor) or with args (Part_Ref)
                elif isinstance(parameter.children[1],
                                (Part_Ref, Structure_Constructor)):
                    field.add_value(key, translate_vertical_dimension(
                            self.vertical_dimension_definition,
                            parameter.children[1].string))

                # For statements with arrays in them (misc_meta_data /
                # synonyms)
                elif isinstance(parameter.children[1], Array_Constructor):
                    for array in walk(parameter.children,
                                      types=Section_Subscript_List):
                        if isinstance(array.children[0],
                                      Char_Literal_Constant):
                            # children0.string = "'foo'"
                            inner_key = array.children[0].string[1:-1]
                        else:
                            # children0.string = "foo"
                            inner_key = array.children[0].string

                        field.add_value(
                                key,
                                (inner_key, array.children[1].string[1:-1])
                                )

                else:

                    field.add_value(key, parameter.children[1].string)

            except Exception as error:
                if field.unique_id:
                    self.LOGGER.warning(
                            "Key: %s on field: %s in file: %s is invalid: %s",
                            key,
                            field.unique_id, file_name, error)
                else:
                    self.LOGGER.warning("Key: %s in file: %s is invalid: %s ",
                                        key, file_name, error)

        return field, valid_field

    @staticmethod
    def extract_multi_line_statement(statement: Level_3_Expr) -> [str, str]:
        """Accepts a fparser Level_3_Expr object as input. This object
        represents multi line statements.
        The function reassembles an arbitrary amount of lines into one
        statement
        :param statement: An fparser Level_3_Expr object
        :return key, value: The key and reassembled value as strings
        """
        value = ""
        key = None
        while isinstance(statement.parent, Level_3_Expr):
            for item in statement.parent.children:

                if isinstance(item, Structure_Constructor_2):
                    if not key:
                        key = item.children[0].string

                    # [1:-1] Removes the quotes from the string
                    value += item.children[1].string[1:-1]

                elif isinstance(item, Char_Literal_Constant):
                    value += item.children[0][1:-1]
            statement = statement.parent

        return key, value


def read_enum(path: str):
    """Reads enumerated values from a file. File should contain only one ENUM
    :param path: Path to the file containing the ENUM
    :return enumerated_values: A list of all enumerated values found
    """
    enumerated_values = []
    reader = FortranFileReader(path)
    parse_tree = F2003_PARSER(reader)

    for enum in walk(parse_tree, types=Enumerator_Def_Stmt):
        for item in enum.children[1].children:
            enumerated_values.append(item.children[0].string)
    return enumerated_values
