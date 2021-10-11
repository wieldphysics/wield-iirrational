#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
""" Setup and build the wavestate.collection package, which acts as the metapackage

Packaging guidance may be found at https://packaging.python.org/tutorials/packaging-projects/
"""
from setuptools import setup

# Settings are primarily in setup.cfg
setup(
    # this ensures setuptools is new enough to use setup.cfg
    setup_requires=[
        'setuptools >= 45.0.0',
    ],
)
