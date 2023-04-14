#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import os.path as path
import subprocess


def pandoc_pdf(out_name):
    def exporter(md_name, folder):
        subprocess.call(["pandoc", "-o", out_name, md_name], cwd=folder)

    return exporter


def npm_markdown_pdf(out_name):
    def exporter(md_name, folder):
        subprocess.call(["markdown-pdf", md_name, "-o", out_name], cwd=folder)

    return exporter
