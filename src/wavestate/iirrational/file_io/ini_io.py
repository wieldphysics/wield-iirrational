#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""



def ini_load(fname):
    import configparser
    dummy_section = 'dummy_section'
    with open(fname) as f:
        file_content = '[{}]\n'.format(dummy_section) + f.read()

    cp = configparser.ConfigParser(interpolation = None)
    cp.read_string(file_content)
    base = dict()
    for op in cp.options(dummy_section):
        base[op] = cp.get(dummy_section, op)
    sections = dict()
    for sec in cp.sections():
        d = dict()
        sections[sec] = d
        for op in cp.options(sec):
            d[op] = cp.get(sec, op)
    if set(base.keys()) ^ set(sections.keys()):
        raise RuntimeError((
            "ini file {} has setting which conflicts with section name").format(fname)
        )
    base.update(sections)
    return base


