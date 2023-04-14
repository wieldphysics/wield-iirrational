#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
""" Version checker to run before building a distribution.

This tool checks for consistency between the git repository tag, the _version.py file
in the module, and the version defined in setup.cfg.

This tool expects to be one directory up from the setup.cfg.
"""
import os
import sys
import subprocess
import configparser
import importlib
import re
from pkg_resources import parse_version
from packaging.version import LegacyVersion

fpath = os.path.dirname(os.path.abspath(__file__))
config = configparser.ConfigParser()

# look for setup.cfg one directory up
config.read(os.path.join(fpath, "..", "setup.cfg"))
version_string = config["metadata"]["version"]
modfile = config["tools.check_versions"]["version_file"]

# import the specified file as the 'version' module
spec = importlib.util.spec_from_file_location(
    "version", os.path.join(fpath, "..", modfile)
)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def check_versions():
    exit = None

    # use the version_info as the base, as it is unaffected by the setuptools_scm
    module_version = ".".join(str(v) for v in module.version_info)

    if module_version != version_string:
        print(
            "WARNING: Stated module version different than setup.py version",
            file=sys.stderr,
        )
        print(
            "         '{0}' != '{1}'".format(module.__version__, version_string),
            file=sys.stderr,
        )
        print("Fix version.py and setup.py for consistency")
        exit = -1

    try:
        git_tag = subprocess.check_output(["git", "describe", "--tags", "--always"])
        git_tag = git_tag.strip()
        git_tag = str(git_tag.decode("utf-8"))
        print("version string:", version_string, "  --------  git tag:", git_tag)

        # remove things like release- or devel- or version- from the start
        match_digit = re.search(r"\d", git_tag)
        git_tag = git_tag[match_digit.start() :]
    except subprocess.CalledProcessError:
        pass
    else:
        if "dev" in version_string:
            p_git_tag = parse_version(git_tag)
            # if the version cannot be parsed, then don't check it
            if isinstance(p_git_tag, LegacyVersion):
                return exit
            elif parse_version(git_tag) > parse_version(version_string):
                print(
                    "WARNING: latex git-tag {} is newer than dev version {}".format(
                        git_tag, version_string
                    ),
                    file=sys.stderr,
                )
                exit = -3
        elif not git_tag.endswith(version_string):
            if version_string in git_tag:
                print(
                    "WARNING: latex git-tag has commits since versioning",
                    file=sys.stderr,
                )
                print(
                    "         '{0}' ---> '{1}'".format(version_string, git_tag),
                    file=sys.stderr,
                )
            else:
                print(
                    "WARNING: latex git-tag different than setup.py version",
                    file=sys.stderr,
                )
                print(
                    "         '{0}' != '{1}'".format(version_string, git_tag),
                    file=sys.stderr,
                )
            print(
                "         Perhaps update versions in setup.py and module.version, git commit, then git tag"
            )
            print(
                "         otherwise fix tag if not yet git pushed to remote (see DISTRIBUTION-README.md)"
            )
            exit = -2

    return exit


if __name__ == "__main__":
    exit = check_versions()
    if exit is not None:
        sys.exit(exit)
