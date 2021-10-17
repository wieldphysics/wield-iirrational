#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np


def squeezerec(key, obj):
    if isinstance(obj, dict):
        for k, v in list(obj.items()):
            obj[k] = squeezerec(k, v)
        return obj
    elif isinstance(obj, np.ndarray):
        if obj.dtype.names is None:
            if not obj.shape:
                return squeezerec(None, obj.item())
            if obj.shape == (0,) or obj.shape == (0, 0):
                return None
            else:
                if obj.dtype == object:
                    lst = []
                    for idx, v in enumerate(obj):
                        lst.append(squeezerec(idx, v))
                    return np.array(lst)
                else:
                    return obj.squeeze()
        else:
            d = dict()
            for name in obj.dtype.names:
                d[name] = squeezerec(name, obj[name])
            return d
    else:
        return obj


def load_matlab(fname):
    import scipy.io
    d = dict()
    d = scipy.io.loadmat(
        fname,
        mdict = d,
        squeeze_me = True,
        chars_as_strings = True,
        #mat_dtype = True,
    )
    d = squeezerec(None, d)
    d.pop('__globals__', None)
    d.pop('__header__', None)
    d.pop('__version__', None)
    return d


def desqueezerec(key, obj):
    if obj is None:
        return np.array([])
    elif isinstance(obj, dict):
        #convert to recarray?
        for k, v in list(obj.items()):
            obj[k] = desqueezerec(k, v)
        return obj
    else:
        return obj


def write_matlab(fname, obj):
    import scipy.io
    d = scipy.io.savemat(
        fname,
        mdict = desqueezerec(None, obj),
    )
    return d
