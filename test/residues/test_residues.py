# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from IIRrational.TFmath.residues import (
    ZPK2residues,
    ZPK2residues_scipy,
    residues2ZPK,
)

import numpy as np


@pytest.mark.parametrize(
    "ZPK",
    [
        (((), (-1,), 10)),
        (((), (-1, -2), 10)),
        (
            (
                (),
                (
                    -1,
                    -1,
                ),
                10,
            )
        ),
        (((), (-1, -1, -2), 10)),
        (((), (-1, -1, -1), 10)),
        (((), (-1, -1, -1, -2), 10)),
        (((), (-1, -1, -1, -1), 10)),
        (((), (-1, -1, -1, -1, -2), 10)),
        (((-3,) * 1, (-1,), 10)),
        (((-3,) * 1, (-1, -2), 10)),
        (
            (
                (-3,) * 1,
                (
                    -1,
                    -1,
                ),
                10,
            )
        ),
        (((-3,) * 1, (-1, -1, -2), 10)),
        (((-3,) * 1, (-1, -1, -1), 10)),
        (((-3,) * 1, (-1, -1, -1, -2), 10)),
        (((-3,) * 1, (-1, -1, -1, -1), 10)),
        (((-3,) * 1, (-1, -1, -1, -1, -2), 10)),
        (((-3,) * 2, (-1,), 10)),
        (((-3,) * 2, (-1, -2), 10)),
        (
            (
                (-3,) * 2,
                (
                    -1,
                    -1,
                ),
                10,
            )
        ),
        (((-3,) * 2, (-1, -1, -2), 10)),
        (((-3,) * 2, (-1, -1, -1), 10)),
        (((-3,) * 2, (-1, -1, -1, -2), 10)),
        (((-3,) * 2, (-1, -1, -1, -1), 10)),
        (((-3,) * 2, (-1, -1, -1, -1, -2), 10)),
        (((-3,) * 3, (-1, -2), 10)),
        (((-3,) * 4, (-1, -2), 10)),
        (((-3,) * 5, (-1, -2), 10)),
        (((-3,) * 6, (-1, -2), 10)),
        (((-3,) * 7, (-1, -2), 10)),
        (((-3,) * 8, (-1, -2), 10)),
    ],
)
def test_residues(ZPK):
    print("ZPK: ", ZPK)
    r_sci, zpk_sci = ZPK2residues_scipy(ZPK)
    r_std, zpk_std = ZPK2residues(ZPK)
    print("r_sci: ", r_sci)
    print("r_std: ", r_std)
    print()
    print("zpk_sci: ", zpk_std)
    print("zpk_std: ", zpk_std)

    ZPK_recomp = residues2ZPK(r_std, zpk_std)
    print("ZPK recomp: ", ZPK_recomp)
