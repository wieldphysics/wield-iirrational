#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


def ZPKrep_DF2_noise(
    ZPK,
    RMSin=1,
    ZPKin=[
        (),
        (
            -0.1,
            -0.1,
        ),
        1,
    ],
    N_quant=3e-16,
    F_nyquist=None,
):
    """
    From DCC G0900928 (M. Evans noise report)
    """
    return


def ZPKrep_BQ_noise(
    ZPK,
    RMSin=1,
    ZPKin=[
        (),
        (
            -0.1,
            -0.1,
        ),
        1,
    ],
    N_quant=3e-16,
    F_nyquist=None,
):
    """
    From DCC G0900928 (M. Evans noise report)
    """
    8 * N_quant * np.max(RMS_in, RMS_out)
    return
