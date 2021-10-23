#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Requires pytest to import
"""


import time
from os import path
import os
import pytest
from os import path
import os
from shutil import rmtree
import contextlib
import matplotlib.pyplot as plt

from . import pytest_utilities


def relfile(_file_, *args, fname=None):
    fpath = path.split(_file_)[0]
    post = path.join(*args)
    fpath = path.join(fpath, post)
    # os.makedirs(fpath, exist_ok = True)
    # os.utime(fpath, None)

    if fname is None:
        return fpath
    else:
        return path.join(fpath, fname)


def relfile_test(_file_, request, pre=None, post=None, fname=None):
    """
    Generates a folder specific to pt.test function
    (provided by using the "request" fixture in the test's arguments)
    """
    if isinstance(pre, (list, tuple)):
        pre = path.join(pre)

    testname = request.node.name
    if pre is not None:
        testname = path.join(pre, testname)

    if isinstance(post, (list, tuple)):
        post = path.join(post)
    if post is not None:
        return relfile(_file_, testname, post, fname=fname)
    else:
        return relfile(_file_, testname, fname=fname)


class Timer(object):
    def __init__(self, N):
        self.N = N

    def __iter__(self):
        return iter(range(self.N))

    def __call__(self):
        return self.interval / self.N

    def __float__(self):
        return self.interval / self.N

    def __str__(self):
        time = self()
        if time > 10:
            return "{:.1f}s".format(time)
        elif time > 1:
            return "{:.2f}s".format(time)
        elif time > 0.001:
            return "{:.1f}ms".format(time * 1e3)
        elif time > 1e-6:
            return "{:.1f}us".format(time * 1e6)
        else:
            return "{:.1f}ns".format(time * 1e9)

    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start


def pytest_addoption(parser):
    parser.addoption(
        "--plot",
        action="store_true",
        dest="plot",
        help="Have tests update plots (it is slow)",
    )

    parser.addoption(
        "--do-benchmarks", action="store_true", help="run slow benchmarking tests"
    )

    parser.addoption(
        "--do-stresstest", action="store_true", help="Run slow repeated stress tests"
    )
    parser.addoption(
        "--browser",
        action="store_true",
        dest="browser",
        help="Make full outputs for browser usage",
    )

    parser.addoption(
        "--plot-verbosity",
        action="store",
        type=int,
        default=None,
        dest="plot_verbosity",
        help="Set the verbosity limit for the markdown output",
    )

    parser.addoption(
        "--plotsections",
        action="store",
        type=str,
        default=None,
        dest="plotsections",
        help="Set the regex to search for plotsection names to limit output",
    )

    parser.addoption(
        "--no-preclear",
        action="store_true",
        default=False,
        dest="no_preclear",
        help="Do not preclear tpaths",
    )

    # parser.addoption("--do-benchmarks", action="store_true",
    #    help="run slow benchmarking tests")


@pytest.fixture
def plot(request):
    return request.config.getvalue("--plot")
    return request.config.option.plot


def tpath_raw_make(request):
    if isinstance(request.node, pytest.Function):
        return relfile_test(
            request.node.function.__code__.co_filename, request, "tresults"
        )
    raise RuntimeError("TPath currently only works for functions")


@pytest.fixture
def tpath_preclear(request):
    tpath_raw = tpath_raw_make(request)
    no_preclear = request.config.getvalue("--no-preclear")
    if not no_preclear:
        rmtree(tpath_raw, ignore_errors=True)


@pytest.fixture
def tpath(request):
    tpath_raw = tpath_raw_make(request)

    os.makedirs(tpath_raw, exist_ok=True)
    os.utime(tpath_raw, None)
    return tpath_raw


@pytest.fixture
def closefigs():
    yield
    plt.close("all")


@pytest.fixture
def tpath_join(request):
    tpath_raw = tpath_raw_make(request)
    first_call = True

    def tpathJ(subpath):
        nonlocal first_call
        if first_call:
            os.makedirs(tpath_raw, exist_ok=True)
            os.utime(tpath_raw, None)
            first_call = False
        return path.join(tpath_raw, subpath)

    return tpathJ


@pytest.fixture
def test_trigger():
    run_store = []

    @contextlib.contextmanager
    def fail(call, **kwargs):
        run_store.append(call)

        def call(did_fail):
            do_call = did_fail
            for k, v in kwargs.items():
                if v:
                    do_call = True
                    break

            if do_call:
                for call in run_store:
                    call(fail=did_fail, **kwargs)
                run_store.clear()

        try:
            yield
        except AssertionError:
            call(True)
            raise
        else:
            call(False)

        return

    return fail


@pytest.fixture()
def ic():
    """
    Fixture to provide icecream imports without requiring that the package exist
    """
    try:
        from icecream import ic

        return ic
    except ImportError:
        pass
    try:
        from IPython.lib.pretty import pprint

        return pprint
    except ImportError:
        from pprint import pprint

        return pprint


@pytest.fixture
def browser(request):
    return request.config.getvalue("--browser")


@pytest.fixture
def plotsections(request):
    return request.config.getvalue("plotsections")


@pytest.fixture
def plot_verbosity(request):
    return request.config.getvalue("plot_verbosity")


@pytest.fixture
def pprint(request, tpath_join):
    fname = tpath_join("output.txt")

    # pushes past the dot
    print("---------------:{}:--------------".format(request.node.name))
    with open(fname, "w") as F:

        def pprint(*args, **kwargs):
            pytest_utilities.pprint(*args, F=F, **kwargs)

        yield pprint
