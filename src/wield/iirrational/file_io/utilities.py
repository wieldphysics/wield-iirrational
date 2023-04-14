#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

from wield import declarative
import collections

NOARG = declarative.unique_generator


def subkey_search(fdict, subkey, default=NOARG):
    if subkey is None:
        return fdict
    subdict = fdict
    skeys = subkey.split(".")

    while skeys:
        if not isinstance(subdict, collections.Mapping):
            raise TypeError("Intermediate type not a dictionary")

        try:
            subdict = subdict[".".join(skeys)]
        except KeyError:
            pass
        else:
            break

        for idx in range(1, len(skeys)):
            semikey = ".".join(skeys[:-idx])
            try:
                subdict = subdict[semikey]
            except KeyError:
                pass
            else:
                skeys = skeys[-idx:]
                break
        else:
            if default is NOARG:
                # this is the for-else syntax, only triggered if the for loop never
                # finds a match
                raise KeyError("Could not recursively find {}".format(subkey))
            else:
                return default
    return subdict
