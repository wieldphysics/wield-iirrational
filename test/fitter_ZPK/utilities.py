# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np

from wavestate.iirrational.fitters_ZPK.MRF import check_jac
from wavestate.iirrational import testing


def check_residuals(fitterMRF):
    # axB = mplfigB()
    for idx in reversed(range(len(fitterMRF.parameters))):
        coding = fitterMRF.coding_maps.p2c[idx]
        # this indicates that the derivative is virtualized
        if coding.derivative_deadzoned:
            continue
        # print(coding, idx)
        try:
            val, select = check_jac(fitterMRF, idx, scaling=0.01)
            if np.any(~np.isfinite(val[select])):
                print(idx, fitterMRF.residuals_jacobian[idx])
                print(idx, val)
            testing.assert_almost_equal(val[select], 1, 2)
        except AssertionError:
            try:
                val, select = check_jac(fitterMRF, idx, scaling=0.0001)
                testing.assert_almost_equal(val[select], 1, 2)
            except AssertionError:
                coding = fitterMRF.coding_maps.p2c[idx]
                print(
                    idx,
                    coding,
                    coding.roots(),
                    coding.reduce(),
                )
                raise
        # axB.ax0.plot(fitterMRF.F_Hz, abs(val))
