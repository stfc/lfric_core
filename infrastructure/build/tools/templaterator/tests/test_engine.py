##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Tests Templaterator functionality.
"""
from hashlib import md5
from pathlib import Path
from textwrap import dedent
from typing import Dict, Optional

from .. import engine


def test_main(tmp_path: Path):
    """
    Ensures application functionality.
    """
    template_file = tmp_path / 'some.f90.template'
    template_text = dedent("""
        module test_{{label}}_mod
          {{type}}({{kind}}) :: variable
        end module test_{{label}}_mod
    """)
    template_file.write_text(template_text)
    keyed_values: Dict[str, Optional[str]] = {
        'label': 'wibble',
        'type': 'real',
        'kind': 'real64'
    }
    output_pattern = str(tmp_path / 'some_{{kind}}.f90')
    engine.main(template_file, keyed_values, output_pattern)

    # Ensure input file is unmodified
    template_hash = md5(template_text.encode()).hexdigest()
    reread_hash = md5(template_file.read_text().encode()).hexdigest()
    assert template_hash == reread_hash

    # Ensure generated file is correct
    generated_file = tmp_path / 'some_real64.f90'
    assert generated_file.exists()
    assert dedent("""
        module test_wibble_mod
          real(real64) :: variable
        end module test_wibble_mod
    """) == generated_file.read_text()
