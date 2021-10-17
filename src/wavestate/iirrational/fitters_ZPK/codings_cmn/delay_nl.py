# -*- coding: utf-8 -*-
"""
"""


import numpy as np
from .base import (
    CodingType,
    I2pi
)


class CodingDelayNL(CodingType):
    N_parameters = 1
    p_delay_nl = 0
    delay_s_min = 0
    delay_s_diff = 0
    deriv_deadzone = 1e-2

    def __init__(self, sys):
        super(CodingDelayNL, self).__init__(sys)

    def clone(self, sys):
        cpy = super(CodingDelayNL, self).clone(sys)
        cpy.delay_s = self.delay_s
        return cpy

    def setup(self, delay_s = None):
        if delay_s is not None:
            self.delay_s = delay_s
        return

    def update(self, delay_nl = None):
        if self.N_parameters > 0:
            self.p_delay_nl = delay_nl
        else:
            assert(delay_nl is None)

    def reduce(self):
        if self.N_parameters > 0:
            return [self.p_delay_nl]
        else:
            return []

    @property
    def N_parameters(self):
        self.delay_s_min = self.sys.delay_s_min
        self.delay_s_diff = self.sys.delay_s_max - self.sys.delay_s_min

        if self.delay_s_diff <= 0:
            N_parameters = 0
            self.delay_s_diff = 0
        else:
            N_parameters = 1
        return N_parameters

    @property
    def delay_s(self):
        return self.delay_s_min + self.delay_s_diff * (1 - np.cos(self.p_delay_nl))/2

    @delay_s.setter
    def delay_s(self, val):
        if self.delay_s_diff <= 0:
            self.p_delay_nl = 0
            return
        if val > self.delay_s_min:
            diffrat = (val - self.delay_s_min) / self.delay_s_diff
            if diffrat > 1:
                self.p_delay_nl = np.pi
            else:
                self.p_delay_nl = np.arccos(1 - diffrat * 2)
        else:
            self.p_delay_nl = 0

    def transfer(self):
        return np.exp(-I2pi * self.delay_s * self.sys.F_Hz)

    @property
    def derivative_deadzoned(self):
        sD = np.sin(self.p_delay_nl)
        if abs(sD) < self.deriv_deadzone:
            return True
        return False

    def derivative(self):
        if self.N_parameters == 0:
            return []
        sD = np.sin(self.p_delay_nl)
        if sD > 0:
            if sD < self.deriv_deadzone:
                sD = self.deriv_deadzone
        else:
            if sD > -self.deriv_deadzone:
                sD = -self.deriv_deadzone
        dDdP = self.delay_s_diff * sD / 2
        return [-I2pi * self.sys.F_Hz * dDdP * np.exp(-I2pi * self.delay_s * self.sys.F_Hz)]

    def roots(self):
        return []

