# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from wavestate.iirrational.testing import IIRrational_data
from wavestate.iirrational import testing

from wavestate.iirrational.fitters_ZPK import ZPKrep2MRF
from wavestate.iirrational.fitters_ZPK import codings_s

from utilities import check_residuals


@pytest.mark.parametrize("coding_map", sorted(codings_s.coding_maps.keys()))
def test_NL_fitter1(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    data = IIRrational_data("testZPK4x")
    print(data.rep_s.zeros)
    data.rep_s.zeros = 0.99 * data.rep_s.zeros
    data.rep_s.poles = 0.99 * data.rep_s.poles
    print(data.rep_s.zeros)
    fitterMRF = ZPKrep2MRF(
        data.rep_s.ZPKrep,
        coding_map=codings_s.coding_maps[coding_map],
    )
    testing.assert_almost_equal(data.rep_s.xfer_fit / fitterMRF.xfer_fit, 1, 4)

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    # axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)


@pytest.mark.parametrize("coding_map", sorted(codings_s.coding_maps.keys()))
def test_NL_fitter2(coding_map):
    data = IIRrational_data("testZPK4x")
    print("Z:", data.rep_s.zeros)
    print("P:", data.rep_s.poles)
    data.rep_s.zeros = 1.01 * data.rep_s.zeros
    data.rep_s.poles = 0.99 * data.rep_s.poles
    fitterMRF = ZPKrep2MRF(
        data.rep_s.ZPKrep,
        coding_map=codings_s.coding_maps[coding_map],
    )

    testing.assert_almost_equal(data.rep_s.xfer_fit / fitterMRF.xfer_fit, 1, 4)

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    assert fitterMRF.residuals_average < 10

    # axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)


@pytest.mark.parametrize("coding_map", sorted(codings_s.coding_maps.keys()))
def test_NL_fitter1_DL(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    # set number 1 is logarithmic data
    data = IIRrational_data("testZPK4xHQ", set_num=1)
    data.rep_s.zeros = 0.99 * data.rep_s.zeros
    data.rep_s.poles = 0.99 * data.rep_s.poles
    print(data.rep_s.zeros)
    fitterMRF = ZPKrep2MRF(
        data.rep_s.ZPKrep,
        coding_map=codings_s.coding_maps[coding_map],
        distance_limit_auto=1,
    )
    testing.assert_almost_equal(data.rep_s.xfer_fit / fitterMRF.xfer_fit, 1, 4)
    fitterMRF.distance_limit_auto = 1

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    # axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)


@pytest.mark.parametrize("coding_map", sorted(codings_s.coding_maps.keys()))
def test_NL_fitter1_DL2(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    # set number 1 is logarithmic data
    data = IIRrational_data("testZPK4xHQ", set_num=1)
    data.rep_s.zeros = 0.99 * data.rep_s.zeros
    data.rep_s.poles = 0.99 * data.rep_s.poles
    print(data.rep_s.zeros)
    fitterMRF = ZPKrep2MRF(
        data.rep_s.ZPKrep,
        coding_map=codings_s.coding_maps[coding_map],
        distance_limit_auto=1,
    )
    testing.assert_almost_equal(data.rep_s.xfer_fit / fitterMRF.xfer_fit, 1, 4)
    fitterMRF.distance_limit_auto = 2

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    # axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)
