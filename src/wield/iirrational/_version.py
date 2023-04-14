#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

version_info = (0, 0, 3, "dev0")
version = ".".join(str(v) for v in version_info)
__version__ = version

# this section of code is removed in "release" branch versions
__version__ = __version__ + "-<unknown>"

try:
    import setuptools_scm
    from pkg_resources import parse_version

    scm_version = setuptools_scm.get_version(
        relative_to=__file__,
        root="../../../",
        fallback_version=__version__,
        version_scheme="guess-next-dev",
    )

    scm_v = parse_version(parse_version(scm_version).base_version)
    version_v = parse_version(parse_version(version).base_version)

    if scm_v != version_v:
        import warnings

        warnings.warn(
            "git base version {} is different than the stored version {}".format(
                scm_v, version_v
            )
        )

    version = scm_version
    __version__ = scm_version
except ValueError as e:
    import warnings

    warnings.warn(str(e))
except (ModuleNotFoundError, TypeError, LookupError):
    pass
