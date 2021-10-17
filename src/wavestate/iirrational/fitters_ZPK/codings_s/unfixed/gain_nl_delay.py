# -*- coding: utf-8 -*-
"""
"""


import numpy as np
from .base import (
    CodingType,
    Ipi,
    I2pi
)


class CodingGainNLDelay(CodingType):
    N_parameters      = 2
    p_gain            = 1
    lock_gain         = False
    hide_gain         = False
    lock_delay        = False
    hide_delay        = False
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

    def setup(
            self,
            hide_all   = None,
            lock_gain  = None,
            hide_gain  = None,
            lock_delay = None,
            hide_delay = None,
            gain       = None,
            delay_s    = None,
    ):
        if hide_all:
            hide_gain = True
            hide_delay = True

        N_parameters = 2
        if lock_gain is not None:
            self.lock_gain = lock_gain
        if hide_gain is not None:
            self.hide_gain = hide_gain
        if self.hide_gain:
            N_parameters -= 1

        if lock_delay is not None:
            self.lock_delay = lock_delay
        if hide_delay is not None:
            self.hide_delay = hide_delay
        if self.hide_delay:
            N_parameters -= 1

        self.N_parameters = N_parameters

        if gain is not None:
            self.p_gain = gain

        if delay_s is not None:
            #self.p_delay = delay_s
            #doesn't actually know the nyqyist, so cant get samples, yet
            #__X = (0 - self.min_delay_samples) / ((self.max_delay_samples - self.min_delay_samples)/2)
            #p_nl_delay = -((__X-1)**2 / (1 - (__X - 1)**2))**.5
            pass

        return

    def update(self, A = None, B = None):
        if self.disable:
            assert(A is None and B is None)
            return
        if self.hide_delay:
            assert(B is None)
            if self.hide_gain:
                assert(A is None)
            else:
                gain = A
                if not self.lock_gain and gain is not None:
                    self.p_gain = gain
        else:
            if self.hide_gain:
                delay = A
                if not self.lock_delay and delay is not None:
                    self.p_nl_delay = delay
            else:
                gain, delay = A, B
                if not self.lock_gain and gain is not None:
                    self.p_gain = gain

                if not self.lock_delay and delay is not None:
                    self.p_nl_delay = delay

    def reduce(self):
        if self.disable:
            return []
        if not self.hide_gain:
            if not self.hide_delay:
                return [
                    self.p_gain,
                    self.p_nl_delay
                ]
            else:
                return [self.p_gain]
        else:
            if not self.hide_delay:
                return [
                    self.p_nl_delay
                ]
            else:
                return []

    @property
    def delay_samples(self):
        delay_samples = self.min_delay_samples + (1 + self.p_nl_delay / (1. + self.p_nl_delay**2)**.5) * ((self.max_delay_samples - self.min_delay_samples)/2)
        return delay_samples

    @property
    def delay_s(self):
        self.delay_samples * self.F_nyquist_Hz

    def delay_set(self, sys, delay_s):
        self.F_nyquist_Hz = sys.F_nyquist_Hz
        assert(False)

    def transfer(self):
        self.F_nyquist_Hz = sys.F_nyquist_Hz
        #gain term (can only be 1 and at the start)
        return self.p_gain * np.exp(-Ipi * self.delay_samples / sys.F_nyquist_Hz * sys.F_Hz)

    def derivative(self):
        self.F_nyquist_Hz = sys.F_nyquist_Hz
        if self.disable:
            return []
        sqp5 = (1. + self.p_nl_delay**2 )**.5
        dDelaydNL = ((self.max_delay_samples - self.min_delay_samples)/2) * (1 / sqp5 - self.p_nl_delay**2 / (1. + self.p_nl_delay**2)**1.5)
        if not self.hide_gain:
            if not self.hide_delay:
                return [
                    np.exp(-Ipi * self.delay_samples / sys.F_nyquist_Hz * sys.F_Hz) if not self.lock_gain else 0,
                    dDelaydNL * self.p_gain * -Ipi * sys.F_Hz / sys.F_nyquist_Hz * np.exp(-Ipi * self.delay_samples * sys.F_Hz / sys.F_nyquist_Hz) if not self.lock_delay else 0,
                ]
            else:
                return [
                    np.exp(-Ipi * self.delay_samples / sys.F_nyquist_Hz * sys.F_Hz) if not self.lock_gain else 0,
                ]
        else:
            if not self.hide_delay:
                return [
                    dDelaydNL * self.p_gain * -Ipi * sys.F_Hz / sys.F_nyquist_Hz * np.exp(-Ipi * self.delay_samples * sys.F_Hz / sys.F_nyquist_Hz) if not self.lock_delay else 0,
                ]
            else:
                return []

    def roots(self):
        return []


