##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
"""
Command line helper functions.
"""
from typing import Dict, List, Optional, Tuple


def parse_kv(kv_list: List[str]) -> Dict[str, Optional[str]]:
    """
    Create dictionary of key, value pairs from argparse list.

    Args:
        kv_list: List of key value pair strings separated by '='.
    return:
        Dictionary of key/value pairs.
    """
    results: Dict[str, Optional[str]] = {}
    if kv_list:
        for item in kv_list:
            key, value = set_kv(item)
            if key in results:
                raise Exception(f"Duplicate key: {key}")
            results[key] = value
    return results


def set_kv(item: str) -> Tuple[str, Optional[str]]:
    """
    Parse the key, value pair.

    Args:
        item: String of key and value separated by '='.
    return:
        Key and value from string
    """
    pair = item.split('=')
    key = pair[0].strip()
    value = None
    if len(pair) > 1:
        value = '='.join(pair[1:])
    return key, value
