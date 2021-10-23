#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from numpy.distutils.system_info import get_info, NotFoundError

from numpy.distutils.misc_util import Configuration


def configuration(parent_package="", top_path=None):
    lapack_opt = get_info("lapack_opt")
    config = Configuration("lapack", parent_package=parent_package, top_path=top_path)

    if not lapack_opt:
        raise NotFoundError("no lapack/blas resources found")

    config.add_extension(
        name="lapack_svd",
        sources=["lapack_svd.pyf.src"],
        extra_info=lapack_opt,
    )
    return config


if __name__ == "__main__":
    from setuptools import setup
    from numpy.distutils.core import setup

    setup(
        description="LAPACK Wrappers",
        author="Lee McCuller",
        author_email="Lee.McCuller@gmail.com",
        **configuration(top_path="").todict()
    )
