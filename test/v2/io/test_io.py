# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest
import os
import numpy as np
import os.path as path

from wavestate.iirrational.testing import IIRrational_data
from wavestate.iirrational.utilities.print import pprint
from wavestate.iirrational.v2.__main__ import main as v2_main

try:
    from IIRrational_test_data import matlab as matlab_test_data
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False

localpath = path.split(__file__)[0]
skipdeco = pytest.mark.skipif(
    module_import_skip, reason="cannot import IIRrational_test_data"
)


def gen_outname(base):
    fname_out = path.join(localpath, base)
    try:
        os.unlink(fname_out)
    except OSError:
        pass
    return fname_out


@skipdeco
def test_mat_datamissing():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    fname_out = path.join(localpath, "out.mat")
    with pytest.raises(SystemExit):
        v2_main(["-D", "test,D", "--mode", "fit", fname, fname_out])


@skipdeco
def test_mat_help():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    with pytest.raises(SystemExit):
        v2_main(["-h"])


@skipdeco
def test_mat_printdata():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    fname_out = path.join(localpath, "out.mat")
    with pytest.raises(SystemExit):
        v2_main(["-D", "data,D", "--mode", "printdata", fname, fname_out])


# TODO test printsettings
# TODO test printconfig


@skipdeco
def test_mat_fit():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    fname_out = gen_outname("out.json")
    v2_main(
        [
            "-D",
            "data,D",
            "--mode",
            "fit",
            fname,
            fname_out,
            "--choose",
            "baseline",
            "-p=-2",
        ]
    )

    fname_out2 = gen_outname("out2.yaml")
    v2_main(
        ["-D", "data,D", "--mode", "fit", fname_out, fname_out2, "--choose", "baseline"]
    )

    fname_out3 = gen_outname("out3.mat")
    v2_main(
        [
            "-D",
            "data,D",
            "--mode",
            "fit",
            fname_out2,
            fname_out3,
            "--choose",
            "baseline",
        ]
    )

    fname_out4 = gen_outname("out4.h5")
    v2_main(
        [
            "-D",
            "data,D",
            "--mode",
            "fit",
            fname_out3,
            fname_out4,
            "--choose",
            "baseline",
        ]
    )

    fname_out5 = gen_outname("out5.pkl")
    v2_main(
        [
            "-D",
            "data,D",
            "--mode",
            "fit",
            fname_out4,
            fname_out5,
            "--choose",
            "baseline",
        ]
    )


@skipdeco
def test_csv_fit():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    from wavestate.iirrational.file_io import load_any

    fdict = load_any(fname=fname, ftype="mat")
    data = fdict["data"]
    darr = np.vstack(
        [
            data["F_Hz"],
            data["xfer"].real,
            data["xfer"].imag,
            data["SNR"],
        ]
    )
    fname_in = path.join(localpath, "in.csv")
    np.savetxt(fname_in, darr, delimiter=",")
    fname_out = gen_outname("out.mat")

    v2_main(
        [
            "-D",
            "data,D",
            "--mode",
            "fit",
            fname_in,
            fname_out,
            "--choose",
            "baseline",
            "-p=-30",
        ]
    )


@skipdeco
def test_mat_full():
    fname = matlab_test_data.matfiles["iir_test_struct.mat"]
    fname_out = gen_outname("out.mat")
    v2_main(
        [
            fname,
            fname_out,
            "--choose",
            "baseline",
            "--relative_degree_min",
            "-2",
        ]
    )
