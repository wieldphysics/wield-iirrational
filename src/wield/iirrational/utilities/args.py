#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

try:
    from collections.abc import Mapping as MappingABC
except ImportError:
    from collections import Mapping as MappingABC


def argscan(*args, **kwargs):
    """
    The attibute name must be given as keyword argument "arg". Giving the name
    allows it to be looked-up in the argument dictionary or locals().
    """
    argname = kwargs.get("arg", None)
    for val in args:
        if val is REQ:
            raise TypeError(
                "Argument {0} is required and must be supplied to above function".format(
                    argname
                )
            )
        elif (argname is not None) and isinstance(val, MappingABC):
            m_arg = val.get(argname, UNSPEC)
            if m_arg is not UNSPEC:
                return m_arg
        elif val is not UNSPEC:
            return val
    # return the last value, likely None or the loop would have returned it
    return val


# Special unique value to specify that arguments are required or not yet specified
REQ = (argscan, "REQUIRED")
UNSPEC = (argscan, "UNSPECIFIED")
