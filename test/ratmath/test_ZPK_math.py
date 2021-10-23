# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from IIRrational.TFmath import ZPKCalc, ratmath
import numpy as np

test_pairs = [
    (
        ((), (-1,), 1),
        ((), (), 1),
    ),
    (
        ((), (-1, -1), 1),
        ((), (), 1),
    ),
    (
        ((), (-1, -1, -1), 1),
        ((), (), 1),
    ),
    (
        ((), (-2, -1, -1), 1),
        ((), (), 1),
    ),
    (
        ((-3,), (-1,), 1),
        ((), (), 1),
    ),
    (
        ((-2, -3), (-1, -1), 1),
        ((), (), 1),
    ),
    (
        ((), (-1, -1, -1), 1),
        ((-3,), (-2,), 1),
    ),
    (
        ((), (-2, -1, -1), 1),
        ((), (), 1),
    ),
    (
        ((), (-2, -1 + 1j, -1 - 1j), 1),
        ((-2 + 1j, -2 - 1j), (-3 + 1j, -3 - 1j), 1),
    ),
]


@pytest.mark.parametrize("ZPK1, ZPK2", test_pairs)
def test_ZPK_sum(ZPK1, ZPK2):
    print("1: ", ZPK1)
    print("2: ", ZPK2)
    ZPK1 = ZPKCalc(*ZPK1, F_nyquist_Hz=None)
    ZPK2 = ZPKCalc(*ZPK2, F_nyquist_Hz=None)
    ZPK3 = ZPKCalc(
        *ratmath.ZPKsum(
            ZPK1,
            ZPK2,
            # scale = 2
        ),
        F_nyquist_Hz=None
    )
    print("3: ", ZPK3)

    F_Hz = np.linspace(0, 100, 100)
    TF1 = ZPK1.TFeval(F_Hz=F_Hz)
    TF2 = ZPK2.TFeval(F_Hz=F_Hz)
    TF3 = ZPK3.TFeval(F_Hz=F_Hz)

    ratio = TF3 / (TF1 + TF2)
    np.testing.assert_almost_equal(ratio, 1)
