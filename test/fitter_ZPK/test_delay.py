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
def test_delay_fitter(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    data = IIRrational_data("delay")
    fitterMRF = ZPKrep2MRF(
        data.rep_s.ZPKrep,
        coding_map=codings_s.coding_maps[coding_map],
    )
    fitterMRF.delay_s = 0.000000
    fitterMRF.codings_revision += fitterMRF.codings_revision + 1
    testing.assert_almost_equal(data.rep_s.xfer_fit / fitterMRF.xfer_fit, 1, 4)
    print("DELAY MRF: ", fitterMRF.delay_s)
    print("DELAY DATA: ", data.rep_s.delay_s)
    print("DELAY MRF MIN: ", fitterMRF.delay_s_min)
    print("DELAY MRF MAX: ", fitterMRF.delay_s_max)

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    assert fitterMRF.residuals_average < 2
    # axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)
