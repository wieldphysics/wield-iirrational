# -*- coding: utf-8 -*-
"""
"""


import numpy as np
from ..codings_cmn import (
    CodingTypeZ,
)


#TODO, fix this class
class CodingGainNLDelay(CodingTypeZ):
    N_parameters      = 2
    max_delay_samples = 100
    min_delay_samples = -1
    p_nl_delay        = -1000

    #delay_samples = min_delay_samples + (1 + p_nl_delay / (1. + p_nl_delay**2)**.5) * ((max_delay_samples - min_delay_samples)/2)
    #(delay_samples - min_delay_samples) / ((max_delay_samples - min_delay_samples)/2) = X = (1 + p_nl_delay / (1. + p_nl_delay**2)**.5)
    #(X-1)**2 = p_nl_delay**2 / (1 + p_nl_delay**2)
    #(X-1)**2 = (1 - (X - 1)**2) * p_nl_delay**2
    #(X-1)**2 / (1 - (X - 1)**2) = p_nl_delay**2
    __X = (0 - min_delay_samples) / ((max_delay_samples - min_delay_samples)/2)
    p_nl_delay = -((__X-1)**2 / (1 - (__X - 1)**2))**.5

    del __X  # delete the temporary so it isn't in the class

    def setup(self, delay_s = None,):
        if delay_s is not None:
            #self.p_delay = delay_s
            #doesn't actually know the nyqyist, so cant get samples, yet
            #__X = (0 - self.min_delay_samples) / ((self.max_delay_samples - self.min_delay_samples)/2)
            #p_nl_delay = -((__X-1)**2 / (1 - (__X - 1)**2))**.5
            pass

        return

    def update(self, nl_delay):
        self.p_nl_delay = nl_delay

    def reduce(self):
        return [self.p_nl_delay]

    @property
    def delay_samples(self):
        delay_samples = self.min_delay_samples + (1 + self.p_nl_delay / (1. + self.p_nl_delay**2)**.5) * ((self.max_delay_samples - self.min_delay_samples)/2)
        return delay_samples

    @property
    def delay_s(self):
        self.delay_samples * self.F_nyquist_Hz

    def delay_set(self, delay_s):
        self.F_nyquist_Hz = self.sys.F_nyquist_Hz
        assert(False)

    def transfer(self):
        self.F_nyquist_Hz = self.sys.F_nyquist_Hz
        #gain term (can only be 1 and at the start)
        return self.p_gain * np.exp(-Ipi * self.delay_samples / self.sys.F_nyquist_Hz * self.sys.F_Hz)

    def derivative(self):
        self.F_nyquist_Hz = self.sys.F_nyquist_Hz
        sqp5 = (1. + self.p_nl_delay**2 )**.5
        dDelaydNL = ((self.max_delay_samples - self.min_delay_samples)/2) * (1 / sqp5 - self.p_nl_delay**2 / (1. + self.p_nl_delay**2)**1.5)
        return [
            dDelaydNL * -Ipi * self.sys.F_Hz / self.sys.F_nyquist_Hz * np.exp(-Ipi * self.delay_samples * self.sys.F_Hz / self.sys.F_nyquist_Hz) if not self.lock_delay else 0,
        ]

    def roots(self):
        return []


