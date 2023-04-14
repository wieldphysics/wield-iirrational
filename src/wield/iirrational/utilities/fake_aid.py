#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


def fitter_order(fitter):
    return max(len(fitter.zeros), len(fitter.poles))


class FitAidFake(object):
    def __init__(
        self,
    ):
        return

    def hint(self, *args, **kwargs):
        default = kwargs.get("default")
        return default

    def hint_arg(self, func_arg, *args, **kwargs):
        if func_arg is None:
            return self.hint(*args, **kwargs)
        else:
            return func_arg

    def log(self, *args, **kwargs):
        print(*args)

    def log_info(self, *args, **kwargs):
        print(*args)

    def log_warn(self, *args, **kwargs):
        print(*args)

    def log_alert(self, *args, **kwargs):
        print(*args)

    def log_debug(self, *args, **kwargs):
        print(*args)

    def log_rationale(self, *args, **kwargs):
        print(*args)


def ensure_aid(aid=None, **kwargs):
    if aid is None:
        return FitAidFake()
    return aid
