##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Tests the CLI helpers.
"""
from pytest import mark

from .. import cli


@mark.parametrize(
    'input_string, key, value',
    [
        ('cheese=beef', 'cheese', 'beef'),
        ('spoon = teapot', 'spoon', 'teapot'),
        ('key with space = value also having spaces',
         'key with space', 'value also having spaces'),
        ('No value specified', 'No value specified', None),
        ('simple key=value with an = in it',
         'simple key', 'value with an = in it')
    ]
)
def test_set_kv(input_string, key, value):
    """
    Ensures variations of key/value pair are understood.
    """
    result = cli.set_kv(input_string)
    assert key, value == result


def test_parse_kv():
    """
    Ensures multiple key/value pairs work.
    """
    result = cli.parse_kv(['fred=barney', 'wilma=betty'])
    assert {'fred': 'barney', 'wilma': 'betty'}, result


def test_parse_kv_single():
    """
    Ensures a single key/value pair works.
    """
    result = cli.parse_kv(['foo=bar'])
    assert {'foo': 'bar'} == result


def test_parse_kv_empty_list():
    """
    Ensures an empty list results in an empty dictionary.
    """
    result = cli.parse_kv([])
    assert not result
