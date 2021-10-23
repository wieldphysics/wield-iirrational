#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import collections
import sys
import re
import yaml
from wavestate.declarative.utilities.future_from_2 import unicode


def yaml_load(fname):
    # HACK: fix loading number in scientific notation
    #
    # https://stackoverflow.com/questions/30458977/yaml-loads-5e-6-as-string-and-not-a-number
    #
    # An apparent bug in python-yaml prevents it from regognizing
    # scientific notation as a float.  The following is a modified version
    # of the parser that recognize scientific notation appropriately.
    yaml_loader = yaml.SafeLoader
    yaml_loader.add_implicit_resolver(
        "tag:yaml.org,2002:float",
        re.compile(
            """^(?:
        [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
        |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
        |\\.[0-9_]+(?:[eE][-+][0-9]+)?
        |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
        |[-+]?\\.(?:inf|Inf|INF)
        |\\.(?:nan|NaN|NAN))$""",
            re.X,
        ),
        list("-+0123456789."),
    )

    with open(fname, "r") as F:
        fdict = yaml.load(F, Loader=yaml_loader)
    return fdict


def dict_representer(dumper, data):
    return dumper.represent_dict(data.iteritems())


yaml.SafeDumper.add_representer(collections.OrderedDict, dict_representer)


def yaml_write(fname, fdict):
    if sys.version < (3, 4):
        with open(fname, "w") as F:
            yaml.safe_dump(fdict, F)
    else:
        with open(fname, "w", encoding="utf8") as F:
            yaml.safe_dump(fdict, F)
