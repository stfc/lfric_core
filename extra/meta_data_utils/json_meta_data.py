##############################################################################
# (c) Crown copyright 2020 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""This module contains functionality for writing LFRic meta data to disk as a
JSON object"""

import hashlib
import json
import logging
import os
from typing import Dict
from entities import Section

LOGGER = logging.getLogger(__name__)


def write_json_meta(data: Dict[str, Section], path: str):
    """Converts the data object into a JSON string, creates a checksum and
    saves it to disk
    :param data: The sorted meta data
    :param path: Path to where the file is to be saved"""

    file_name = "LFRic_meta_data.json"
    serialise_data = json.dumps(data, default=lambda o: o.__dict__, indent=4)

    data_for_disk = {"meta_data": json.loads(serialise_data)}
    set_checksum(data_for_disk)

    with open(os.path.join(path, file_name), "w") as file_handle:
        LOGGER.info("Writing %s", file_name)
        json.dump(data_for_disk, file_handle, indent=4)
        file_handle.write(os.linesep)


def set_checksum(meta_data: Dict[str, Dict]):
    """Calculate the checksum for meta_data, then adds the value to
    meta_data under the checksum key. meta_data is modified in place.
    :param meta_data: A Dictionary containing the data to be saved"""

    obj_str = json.dumps(meta_data, sort_keys=True)
    checksum = 'md5: ' + hashlib.md5(obj_str.encode('utf-8')).hexdigest()
    meta_data['checksum'] = checksum
