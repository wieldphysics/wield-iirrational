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
import declarative
import numpy as np
from declarative.utilities.future_from_2 import unicode

from ... import fitters_ZPK
from ...utilities import args
from ... import representations


class ArgumentError(ValueError):
    pass


def mapcheck_bool(aid, aname, val):
    if isinstance(val, (str, unicode)):
        if val.lower() in ['true', 'yes', '1']:
            val = True
        elif val.lower() in ['false', 'no', '0']:
            val = False
        else:
            raise ArgumentError((
                "Argument {} has unrecognized bool specifier {}"
            ).format(aname, val))
    return bool(val)

def float_check(aid, aname, val):
    try:
        return float(val)
    except ValueError:
        raise ArgumentError((
            "argument {}={} must be a float"
        ).format(aname, val))

def mapcheck_float(aid, aname, val):
    try:
        return float(val)
    except ValueError:
        raise ArgumentError((
            "argument {}={} must be a float"
        ).format(aname, val))

def mapcheck_positive_float(aid, aname, val):
    val = float(val)
    if not (val > 0):
        raise ArgumentError((
            "argument {}={} must be positive"
        ).format(aname, val))
    return val

def mapcheck_nonnegative_float(aid, aname, val):
    val = float(val)
    if not (val >= 0):
        raise ArgumentError((
            "argument {}={} must not be negative"
        ).format(aname, val))
    return val

def mapcheck_nonnegative_float_orNone(aid, aname, val):
    if isinstance(val, (str, unicode)):
        if val.lower() in ['none', 'null']:
            val = None
    if val is None:
        return val
    return mapcheck_nonnegative_float(aid, aname, val)

def mapcheck_positive_float_orNone(aid, aname, val):
    if isinstance(val, (str, unicode)):
        if val.lower() in ['none', 'null']:
            val = None
    if val is None:
        return val
    return mapcheck_positive_float(aid, aname, val)


def mapcheck_int(aid, aname, val):
    try:
        val = int(val)
    except ValueError:
        raise ArgumentError((
            "argument {}={} must be an integer"
        ).format(aname, val))

    return val


def mapcheck_positive_int(aid, aname, val):
    try:
        val = int(val)
    except ValueError:
        raise ArgumentError((
            "argument {}={} must be an integer"
        ).format(aname, val))

    #TODO error reporting
    if not (val > 0):
        raise ArgumentError((
            "argument {}={} must be positive"
        ).format(aname, val))
    if not (val - float(val) == 0):
        raise ArgumentError((
            "argument {}={} must be an integer"
        ).format(aname, val))
    return val


def mapcheck_nonnegative_int(aid, aname, val):
    try:
        val = int(val)
    except ValueError:
        raise ArgumentError((
            "argument {}={} must be an integer"
        ).format(aname, val))

    #TODO error reporting
    if not (val >= 0):
        raise ArgumentError((
            "argument {}={} must be non-negative"
        ).format(aname, val))
    if not (val - float(val) == 0):
        raise ArgumentError((
            "argument {}={} must be an integer"
        ).format(aname, val))
    return val


def mapcheck_positive_int_orNone(aid, aname, val):
    if isinstance(val, (str, unicode)):
        if val.lower() in ['none', 'null']:
            val = None
    if val is None:
        return val
    return mapcheck_positive_int(aid, aname, val)


def mapcheck_nonnegative_int_orNone(aid, aname, val):
    if isinstance(val, (str, unicode)):
        if val.lower() in ['none', 'null']:
            val = None
    if val is None:
        return val
    return mapcheck_nonnegative_int(aid, aname, val)


def mapcheck_int_orNone(aid, aname, val):
    if isinstance(val, (str, unicode)):
        if val.lower() in ['none', 'null']:
            val = None
    if val is None:
        return val
    return mapcheck_int(aid, aname, val)


def kw_ZPKrep_build(args, kwargs):
    """
    """
    kw = dict()
    ZPKrep = None
    for arg in args:
        if isinstance(arg, collections.Mapping):
            kw.update(arg)
        elif isinstance(arg, (
            representations.ZPKwData,
            fitters_ZPK.MultiReprFilterBase,
        )):
            if ZPKrep is None:
                ZPKrep = arg.ZPKrep
            else:
                raise ArgumentError(
                    "Only a single ZPKrep may be specified"
                    " in positional arguments"
                )
        else:
            raise ArgumentError(
                "Positional Arguments are not allowed,"
                " except for a single master ZPKrep, or a set of"
                " default keyword argument dictionaries"
            )
    kw.update(kwargs)
    return kw, ZPKrep


def grab_kwarg_hints(aid, kw, kwdesc, kwput = None):
    def eval_hint(hint_name):
        kwmeta = kwdesc[hint_name]

        reqs = kwmeta.get('require_hints', None)
        if reqs is not None:
            for hname in reqs:
                if not aid.hint_has(hname):
                    eval_hint(hname)
                if not aid.hint_has(hname):
                    raise RuntimeError("Hint/Kwarg evaluation out of order. BUG!")

        found = _grab_kwargs(aid, kw, kwmeta, hint_name, kwput = kwput)
        if found:
            aid.hint_setdefault(hint_name, next(iter(found.values())))
        #TODO, make better error message for case when default is not given
        try:
            default = kwmeta['default']
        except KeyError:
            pass
        else:
            if callable(default):
                default = default(aid, hint_name)
            aid.hint_setdefault(hint_name, default)

    for hint_name in kwdesc.keys():
        eval_hint(hint_name)


def grab_kwargs(aid, kw, kwdesc, argname, kwput = None):
    kwmeta = kwdesc[argname]
    found = _grab_kwargs(aid, kw, kwmeta, argname, kwput = kwput)
    if found:
        return next(iter(found.values()))

    #TODO, make better error message for case when default is not given
    default = kwmeta['default']
    if callable(default):
        default = default(aid, argname)
    return default


def _grab_kwargs(aid, kw, kwmeta, argname, kwput = None):
    reqs = kwmeta.get('require_hints', None)
    if reqs is not None:
        for hname in reqs:
            if not aid.hint_has(hname):
                raise RuntimeError("Hint/Kwarg evaluation out of order. BUG!")
    drop_val = kwmeta.get('drop_val', args.UNSPEC)
    drop_vals = list(kwmeta.get('drop_vals', []))
    if drop_val is not args.UNSPEC:
        drop_vals.extend(drop_val)

    pop = kwmeta.get('pop', True)
    mapcheck = kwmeta.get('mapcheck', None)
    found = {}
    #use a list to modify it from the function
    discrepancy = [False]

    def check_find(aname):
        if pop:
            val = kw.pop(aname, args.UNSPEC)
        else:
            val = kw.get(aname, args.UNSPEC)

        if val is args.UNSPEC:
            return False

        if kwput is not None:
            kwput[aname] = val

        if val in drop_vals:
            return False

        if mapcheck is not None:
            val = mapcheck(aid, aname, val)

        if found:
            if np.any(next(iter(found.values())) != val):
                discrepancy[0] = True
        found[aname] = val

        return True

    prefname = kwmeta.get('name', argname)
    check_find(prefname)

    for aname in kwmeta.get('aliases', []):
        check_find(aname)

    for aname in kwmeta.get('aliases_bad', []):
        if check_find(aname):
            aid.log_warn(
                2,
                "Argument '{}' is a deprecated alias."
                " Use '{}' instead".format(aname, prefname),
            )

    if discrepancy[0]:
        kv_strs = []
        for k, v in found.items():
            kv_strs.append("{}: {}".format(k, v))
        raise ArgumentError(
            "Inconsistent values given for aliased arguments\n" + '\n'.join(kv_strs)
        )

    return found


def check_remaining_arguments(kw, kwdict):
    names = set()
    alias_map = {}
    bad_map = {}
    for hname, hdict in kwdict.items():
        name = hdict.get('name', hname)
        names.add(name)
        for alias in hdict.get('aliases', []):
            alias_map[alias] = name
        for alias in hdict.get('aliases_bad', []):
            bad_map[alias] = name

    allnames = set(names)
    allnames.update(alias_map.keys())
    allnames.update(bad_map.keys())
    import difflib
    k_lines = dict()
    for k in kw.keys():
        l = []
        k_lines[k] = l
        matches = difflib.get_close_matches(k, allnames, n = 5)
        print("Matches: ", k, matches)
        mainset = set()
        for match in matches:
            if match in names:
                mainset.add(name)
                l.append("'{}'".format(match))
            elif match in alias_map:
                mainmatch = alias_map[match]
                l.append("'{}', alias for '{}'".format(match, mainmatch))
            else:
                mainmatch = bad_map[match]
                l.append("'{}', deprecated alias for '{}'".format(match, mainmatch))
    stack_lines = []
    for k, lines in sorted(k_lines.items()):
        if len(lines) > 0:
            stack_lines.append("'{}':\n\t".format(k) + "\n\t".join(lines))
        else:
            stack_lines.append("'{}': <no similar arguments recognized>".format(k))
    raise ArgumentError((
        "Unrecognized keyword arguments. Listed with potential matches below,"
        "\nor call with argument 'help'=True for more details\n"
        "{}"
    ).format('\n'.join(stack_lines)))
    return


def cplx_iIjJ(val):
    val = val.replace('i', 'j')
    val = val.replace('I', 'j')
    val = val.replace('J', 'j')
    return complex(val)

def cplx_iIjJ_list(val):
    vals = val.split()
    vals2 = []
    for v in vals:
        vals2.extend(v.split(','))
    vals3 = []
    for v in vals2:
        vals3.extend(v.split(';'))
    return [cplx_iIjJ(v) for v in vals3 if v is not '']

def transfer_kw(kwA, kwB, kwdict, pop = False):
    for hname, hdict in kwdict.items():
        if pop:
            name = hdict.pop('name', hname)
        else:
            name = hdict.get('name', hname)
        normalize = hdict.get('normalize', None)
        val = kwA.get(name, declarative.NOARG)
        if val is not declarative.NOARG:
            if normalize is not None:
                val = normalize(val)
            kwB[name] = val

        for aname in hdict.get('aliases', []):
            if pop:
                val = kwA.pop(aname, declarative.NOARG)
            else:
                val = kwA.get(aname, declarative.NOARG)
            if val is not declarative.NOARG:
                if normalize is not None:
                    val = normalize(val)
                kwB[aname] = val

        for aname in hdict.get('aliases_bad', []):
            if pop:
                val = kwA.pop(aname, declarative.NOARG)
            else:
                val = kwA.get(aname, declarative.NOARG)
            if val is not declarative.NOARG:
                if normalize is not None:
                    val = normalize(val)
                kwB[aname] = val
