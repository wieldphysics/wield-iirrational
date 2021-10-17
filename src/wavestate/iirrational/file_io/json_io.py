# -*- coding: utf-8 -*-
"""
"""

import json
import sys


def load_json(fname):
    with open(fname) as F:
        fdict = json.load(F)
    return fdict

def write_json(fname, fdict):
    if sys.version_info < (3, 4):
        with open(fname, 'w') as F:
            json.dump(
                fdict, F,
                indent=4,
                ensure_ascii = False
            )
    else:
        with open(fname, 'w', encoding='utf8') as F:
            json.dump(
                fdict, F,
                indent=4,
                ensure_ascii = False
            )
    return
