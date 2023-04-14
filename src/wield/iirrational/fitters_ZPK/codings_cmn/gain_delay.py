#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import numpy as np
from .base import (
    CodingType,
    # Ipi,
    I2pi,
)


class CodingGain(CodingType):
    N_parameters = 1
    p_gain = 1

    def __init__(self):
        self.coding_id = 0
        return

    def clone(self, sys):
        new = self.__class__()
        new.update(gain=self.p_gain)
        return new

    def setup(self, gain=None):
        if gain is not None:
            self.p_gain = gain
        return

    @property
    def gain(self):
        return self.p_gain

    @gain.setter
    def gain(self, gain):
        self.p_gain = gain

    def update(self, gain):
        self.p_gain = gain

    def reduce(self):
        return [self.p_gain]

    def transfer(self):
        # gain term (can only be 1 and at the start)
        return self.p_gain

    def derivative(self):
        return [1]

    def roots(self):
        return []


class CodingDelay(CodingType):
    N_parameters = 1
    p_delay_s = 0

    def setup(self, delay_s=None):
        if delay_s is not None:
            self.p_delay_s = delay_s
        return

    @property
    def delay_s(self):
        return self.p_delay_s

    @delay_s.setter
    def delay_s(self, val):
        self.p_delay_s = val

    def update(self, delay):
        self.p_delay_s = delay

    def reduce(self):
        return [self.p_delay_s]

    def transfer(self):
        return np.exp(-I2pi * self.p_delay_s * self.sys.F_Hz)

    def derivative(self):
        return [-I2pi * self.sys.F_Hz * np.exp(-I2pi * self.p_delay_s * self.sys.F_Hz)]

    def roots(self):
        return []
