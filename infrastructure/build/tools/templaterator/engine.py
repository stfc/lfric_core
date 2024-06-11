##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Reifies compilable source from an abstract template of such.

Keyed values are inserted into the template at specified locations to produce
something a compiler can understand.
"""
from pathlib import Path
from typing import Dict, Optional

from jinja2 import BaseLoader, Environment, FileSystemLoader


def main(source_path: Path,
         kv_dict: Dict[str, Optional[str]],
         output_file: str) -> None:
    """
    Main method

    Args:
        source_path: String for path to template source.
        kv_dict: List of dictionaries to match and replace.
        output_file: Pattern for output filename with template.
    """
    environment = Environment(
        variable_start_string='{{', variable_end_string='}}',
        loader=FileSystemLoader(source_path.parent),
        keep_trailing_newline=True
    )
    template = environment.get_template(source_path.name)
    file_name_template = Environment(
        loader=BaseLoader()).from_string(output_file)
    filename = file_name_template.render(kv_dict)
    content = template.render(kv_dict)

    with open(filename, mode="tw", encoding="utf-8") as message:
        message.write(content)
        print(f"... wrote {filename}")
