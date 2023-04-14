#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import re

_re_any_words = re.compile("\s*\S+")
_re_leading_space = re.compile("\s*")


def padding_remove(docstring_like, tabsize=4):
    """
    Both removes the leading line paddings on triple quoted strings and removes common indentation
    (expands tabs to do so)
    """
    docstring_like = docstring_like.strip()
    docstring_like = docstring_like.expandtabs(tabsize)
    docstring_lines = docstring_like.splitlines(True)

    spacings = []
    for line in docstring_lines[1:]:
        if len(line.strip()) == 0:
            continue
        spacings.append(len(_re_leading_space.match(line).group(0)))
    if spacings:
        common_spacing = min(spacings)
    else:
        common_spacing = 0
    common_tab_spacing = (common_spacing // tabsize) * tabsize

    if docstring_lines:
        fixed_lines = [docstring_lines[0]]
        for line in docstring_lines[1:]:
            fixed_lines.append(line[common_tab_spacing:])
    else:
        fixed_lines = []

    return "".join(fixed_lines)
