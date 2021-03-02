##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""This module contains functionality for writing LFRic meta data to disk as a
rose-meta.conf file"""
import logging
import os
from textwrap import wrap
from typing import Dict

from entities import Section

LOGGER = logging.getLogger("lfric_meta_data_parser")


def write_file(path: str, file_name, data: str):
    """
    Writes rose suite files to disk
    :param path: The path where the file is to be written
    :param file_name: The name of the file to be written
    :param data: file as a string
    """
    with open(path + file_name + ".conf", 'w') as file:
        file.write(data)


def create_rose_meta(meta_data: Dict, directory: str, file_name):
    """
    Creates a rose_meta.conf file using the supplied meta data
    :param meta_data: Dict containing all meta data
    :param directory: The directory the file will be saved in
    :param file_name: The name of the file that the meta data will be
    written to
    :return rose_meta: The rose-meta file as a string, ready to be written to
    disk
    """
    LOGGER.info("Creating %s.conf", file_name)
    rose_meta = ""

    rose_meta += """
[field_config]
title=LFRic Field Configuration
"""

    for section in meta_data["sections"].values():
        rose_meta += f"""
[field_config:{section.name}]
title={section.title}"""

        for group in section.groups.values():
            rose_meta += f"""
[field_config:{section.name}:{group.name}]
title={group.title}"""

            for field in group.fields.values():
                rose_meta += f"""
[field_config:{section.name}:{group.name}={field.unique_id}]
type=boolean
title=Enable {field.item_title}
trigger=field_config:{section.name}:{group.name}={field.unique_id}{
                field.trigger}
help=Unit of Measure: {field.units}
    =Function Space: {field.function_space}
    =Data type: {field.data_type}
    =Time step: {field.time_step}
    =Interpolation: {field.recommended_interpolation}
"""
                if field.vertical_dimension:
                    attribute_string = f"    =vertical_dimension:{os.linesep}"
                    for key, value in field.vertical_dimension.items():
                        attribute_string += f"       ={key}: " \
                                            f"{str(value)}{os.linesep}"

                    rose_meta += attribute_string

                if field.synonyms:
                    attribute_string = f"    =Synonyms:{os.linesep}"
                    for key, values in field.synonyms.items():
                        for value in values:
                            attribute_string += f"       ={key.value}: " \
                                                f"{str(value)}{os.linesep}"

                    rose_meta += attribute_string

                rose_meta += f"""
description={wrap(field.description, width=100)}
           =For more information on {field.item_title}, see the help text

[field_config:{section.name}:{group.name}={field.unique_id}__checksum]
type=boolean
title=Enable Checksum for {field.item_title} checksum
compulsory=true
"""

    rose_meta = add_file_meta(meta_data, rose_meta)
    rose_meta = add_vertical_meta(rose_meta,
                                  meta_data["standard_level_markers"])

    write_file(directory + "meta/", file_name, rose_meta)


def add_file_meta(meta_data: Dict, rose_meta: str) -> str:
    """Adds meta data for file output in rose.
    :param meta_data: Dict containing parsed meta data
    :param rose_meta: String that the file output data will be appended to
    :return rose_meta: String with file output data appended to it"""

    values_list = []
    titles_list = []

    for section in meta_data["sections"].values():
        for group in section.groups.values():
            for field in group.fields.values():

                values_list.append(field.unique_id)
                titles_list.append(
                        section.title + ": " + group.title +
                        ": " + field.item_title)

    values = ', '.join(values_list)
    titles = '", "'.join(titles_list)

    rose_meta += f"""[output_stream]
duplicate=true
macro=add_section.AddField, add_section.AddStream

[output_stream=name]
type=character

[output_stream=timestep]
type=character

[output_stream:field]
duplicate=true
macro=add_section.AddField

[output_stream:field=id]
values={values}
value-titles="{titles}"

[output_stream:field=temporal]
values=instant,average,accumulate,minimum,maximum,once"""
    return rose_meta


def add_vertical_meta(rose_meta: str, levels) -> str:
    """Adds data about vertical dimensions. Currently static, will be further
    developed in the future
    :param rose_meta: String that the vertical dimension data will be appended
    to
    :param levels: A List of model levels
    :return rose_meta: String with vertical dimension data appended to it"""
    rose_meta += """
[vertical_dimension]
duplicate=true
title=Vertical Dimension

[vertical_dimension=name]
title=Name
description=Name of the vertical dimension
help=The name used to identify this vertical dimension when associating a field
     with it in Rose
type=character
compulsory=true
fail-if=len(this) == 0 # Name must be specified
sort-key=01

[vertical_dimension=positive]
title=Positive
description=The positive direction
help=The positive direction of the vertical axis, either up or down
values=up, down
compulsory=true
sort-key=02

[vertical_dimension=units]
title=Units
description=Unit of measure
help=The unit of measure for this vertical axis is restricted to be in metres
values=m
compulsory=true
sort-key=03

[vertical_dimension=level_definition]
title=Level boundaries
description=Boundaries of levels in ascending order
help=Positive numbers defining the edges of each level in the vertical
     dimension. The boundaries should be entered in ascending order
length=:
type=real
macro=level_definition.Validator, level_definition.Transformer
range=0:
fail-if=len(this)<2 # There must be at least two level boundaries
compulsory=true
sort-key=04
"""
    num = 1001
    for level in levels:
        rose_meta += f"""[vertical_dimension={level}]
title={level.replace("_", " ".title())}
description=A Model Level
type=integer

range=0:
# Layer out of range
fail-if=this > len(vertical_dimension=level_definition)-1;
sort-key=model-levels-{num}"""
        rose_meta += os.linesep
        num += 1
    return rose_meta
