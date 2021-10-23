#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import time
from os import path
import os
from wavestate.utilities.mpl import mplfigB

try:
    import icecream
except ImportError:
    icecream = None
    pass

try:
    from IPython.lib.pretty import pprint, pretty

    pformat = pretty
except ImportError:
    from pprint import pprint, pformat


def relfile(_file_, *args, fname=None):
    fpath = path.split(_file_)[0]
    post = path.join(*args)
    fpath = path.join(fpath, post)
    os.makedirs(fpath, exist_ok=True)
    os.utime(fpath, None)

    if fname is None:
        return fpath
    else:
        return path.join(fpath, fname)


def relfile_test(_file_, request, pre=None, post=None, fname=None):
    """
    Generates a folder specific to pt.test function
    (provided by using the "request" fixture in the test's arguments)
    """
    if isinstance(pre, (list, tuple)):
        pre = path.join(pre)

    testname = request.node.name
    if pre is not None:
        testname = path.join(pre, testname)

    if isinstance(post, (list, tuple)):
        post = path.join(post)
    if post is not None:
        return relfile(_file_, testname, post, fname=fname)
    else:
        return relfile(_file_, testname, fname=fname)


class Timer(object):
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start


def pprint(*args, F=None, pretty=True, **kwargs):
    outs = []
    if pretty:
        for arg in args:
            outs.append(pformat(arg))
    else:
        outs = args
    if F is not None:
        print(*outs, file=F)
    if icecream is not None:
        icecream.DEFAULT_OUTPUT_FUNCTION(" ".join(outs), **kwargs)
    else:
        print(*outs, **kwargs)
