#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

# import numpy as np

# oldprint = print
# def print(*args):
#    pargs = []
#    for a in args:
#        try:
#            if isinstance(a, np.array):
#                rep = a
#            else:
#                rep = unicode(a, 'utf-8')
#        except TypeError:
#            rep = repr(a)
#        pargs.append(rep)
#    oldprint(*args)

try:
    from IPython.lib.pretty import pprint
except ImportError:
    from pprint import pprint
