# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest
import numpy as np
import os.path as path

from IIRrational.testing import IIRrational_data
from IIRrational import file_io
from IIRrational.file_io import (
    matlab_io,
)

from IIRrational.utilities.print import pprint

try:
    from IIRrational_test_data import matlab as matlab_test_data
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_mat():
    fname = matlab_test_data.matfiles["iir_test_vars.mat"]
    fdict = matlab_io.load_matlab(fname)
    fdict["F_Hz"]
    assert fdict["xfer"].dtype == np.complex
    fdict["SNR"]
    assert fdict["noneA"] is None
    assert fdict["noneC"] is None


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_mat_struct():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    fdict = matlab_io.load_matlab(fname)
    # pprint(fdict['data']['xfer'])
    assert file_io.subkey_search(fdict, "data.xfer").dtype == complex
    # pprint(fdict['data']['confuse'])
    assert fdict["data"]["noneA"] is None
    assert fdict["data"]["noneC"] is None


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_mat_struct_RI():
    fname = matlab_test_data.matfiles["iir_test_struct_RI.mat"]
    fdict = matlab_io.load_matlab(fname)
    # pprint(fdict)


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_mat_struct_RT():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    fdict = matlab_io.load_matlab(fname)
    test_path = path.split(__file__)[0]
    fpath = path.join(test_path, "iir_test_struct_RT.mat")
    matlab_io.save_matlab(fpath, fdict)
    fdict2 = matlab_io.load_matlab(fpath)

    assert np.all(fdict["data"]["F_Hz"] == fdict2["data"]["F_Hz"])
    assert np.all(fdict["data"]["xfer"] == fdict2["data"]["xfer"])
    assert np.all(fdict["data"]["SNR"] == fdict2["data"]["SNR"])
    assert np.all(fdict["data"]["noneA"] == fdict2["data"]["noneA"])
    assert np.all(fdict["data"]["noneC"] == fdict2["data"]["noneC"])
    print(str(fdict["data"]["confuse"]))
    print(str(fdict2["data"]["confuse"]))
    assert str(fdict["data"]["confuse"]) == str(fdict2["data"]["confuse"])
