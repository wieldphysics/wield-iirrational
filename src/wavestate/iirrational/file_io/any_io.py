#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Some notes about writing/reading a common dictionary format across these files:
Some of the issues can only be resolved by code-side normalization functions.

Only dictionaries and arrays of numbers should be supported
 - hdf5 cannot have heterogeneous lists of objects (others can)
 - matlab array squeezing can truncate 1x1 arrays to single elements. Only
   code-side normalization can fix this
 - json/yaml must split complex numbers into dicts with real/imag parts. This
   could/should be reversed on load, but currently requires code normalization
   to reverse
"""

import collections
import numpy as np
import copy

from . import types


def load_any(
    fname,
    ftype,
):
    if ftype == "hdf5":
        from . import hdf5_io

        fdict = hdf5_io.load_hdf5(fname)
    elif ftype == "json":
        from . import json_io

        fdict = json_io.load_json(fname)
    elif ftype == "yaml":
        from . import yaml_io

        fdict = yaml_io.yaml_load(fname)
    elif ftype == "pickle":
        from . import pickle_io

        fdict = pickle_io.load_pickle(fname)
    elif ftype == "ini":
        from . import ini_io

        fdict = ini_io.ini_load(fname)
    elif ftype == "mat":
        from . import matlab_io

        fdict = matlab_io.load_matlab(fname)
    elif ftype == "special":
        try:
            import IIRrational_test_data

            IIRrational_data = IIRrational_test_data.IIRrational_data_dev
        except (ImportError, AttributeError):
            from ..testing import IIRrational_data
        fdict = IIRrational_data(fname)
    else:
        raise RuntimeError("Unsupported filetype")
    return fdict


def cull_None(obj):
    if isinstance(obj, collections.Mapping):
        dels = []
        for k, v in obj.items():
            v2 = cull_None(v)
            if v2 is None:
                dels.append(k)
            else:
                obj[k] = v2
        for k in dels:
            del obj[k]
        return obj
    elif isinstance(obj, np.ndarray):
        if obj.dtype == object:
            for idx, v in np.ndenumerate(obj):
                obj[idx] = cull_None(v)
            return obj
        else:
            return obj
    elif isinstance(obj, list):
        for idx, v in enumerate(obj):
            obj[idx] = cull_None(v)
        return obj
    elif isinstance(obj, tuple):
        obj2 = [None] * len(obj)
        for idx, v in enumerate(obj):
            obj2[idx] = cull_None(v)
        return tuple(obj2)
    return obj


def normalize_ndarray(obj):
    if isinstance(obj, collections.Mapping):
        for k, v in obj.items():
            obj[k] = normalize_ndarray(v)
        return obj
    elif isinstance(obj, np.ndarray):
        return obj
    elif isinstance(obj, np.generic):
        # captures scalars given the previous if
        #return np.asscalar(obj)
        return obj.item()
    elif isinstance(obj, (list, tuple)):
        return np.asarray(obj)
    return obj


def fix_complex(obj):
    normalize_ndarray(obj)
    if isinstance(obj, collections.Mapping):
        for k, v in obj.items():
            obj[k] = fix_complex(v)
        return obj
    elif isinstance(obj, np.ndarray):
        if np.iscomplexobj(obj):
            objD = collections.OrderedDict()
            objD["real"] = obj.real
            objD["imag"] = obj.imag
            objD["__type__"] = "[complex]"
            obj = objD
            return obj
        elif obj.dtype == object:
            for idx, v in np.ndenumerate(obj):
                obj[idx] = fix_complex(v)
            return obj
        else:
            return obj
    elif isinstance(obj, list):
        for idx, v in enumerate(obj):
            obj[idx] = fix_complex(v)
        return obj
    elif isinstance(obj, tuple):
        obj2 = [None] * len(obj)
        for idx, v in enumerate(obj):
            obj2[idx] = fix_complex(v)
        return tuple(obj2)
    elif isinstance(obj, complex):
        objD = collections.OrderedDict()
        objD["real"] = obj.real
        objD["imag"] = obj.imag
        objD["__type__"] = "[complex]"
        obj = objD
    return obj


def fix_ndarray(obj):
    if isinstance(obj, collections.Mapping):
        for k, v in obj.items():
            obj[k] = fix_ndarray(v)
        return obj
    elif isinstance(obj, np.ndarray):
        if obj.dtype == object:
            obj = fix_ndarray(obj.tolist())
            return obj
        else:
            obj = obj.tolist()
            return obj
    elif isinstance(obj, np.generic):
        #return np.asscalar(obj)
        return obj.item()
    elif isinstance(obj, list):
        for idx, v in enumerate(obj):
            obj[idx] = fix_ndarray(v)
        return obj
    elif isinstance(obj, tuple):
        obj2 = [None] * len(obj)
        for idx, v in enumerate(obj):
            obj2[idx] = fix_ndarray(v)
        return tuple(obj2)
    return obj


def write_any(
    fname,
    ftype,
    fdict,
):
    features = types.type2features[ftype]
    fdict_orig = None

    if fdict_orig is None:
        fdict_orig = fdict
        fdict = copy.deepcopy(fdict)
    cull_None(fdict)
    if not features["complex"]:
        if fdict_orig is None:
            fdict_orig = fdict
            fdict = copy.deepcopy(fdict)
        fix_complex(fdict)

    if not features["ndarray"]:
        if fdict_orig is None:
            fdict_orig = fdict
            fdict = copy.deepcopy(fdict)
        fix_ndarray(fdict)

    if ftype == "hdf5":
        from . import hdf5_io

        fdict = hdf5_io.write_hdf5(fname, fdict)
    elif ftype == "mat":
        from . import matlab_io

        fdict = matlab_io.write_matlab(fname, fdict)
    elif ftype == "json":
        from . import json_io

        fdict = json_io.write_json(fname, fdict)
    elif ftype == "yaml":
        from . import yaml_io

        fdict = yaml_io.yaml_write(fname, fdict)
    elif ftype == "pickle":
        from . import pickle_io

        fdict = pickle_io.write_pickle(fname, fdict)
    elif ftype == "ini":
        from . import ini_io

        fdict = ini_io.ini_write(fname, fdict)
    elif ftype == "csv":
        raise RuntimeError("Unsupported Output Type")
    return fdict
