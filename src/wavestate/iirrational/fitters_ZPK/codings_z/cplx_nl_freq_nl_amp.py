# -*- coding: utf-8 -*-
"""
"""


import numpy as np

from ..codings_cmn import (
    CodingTypeZ,
    Ipi,
)


class CodingnlFnlA(CodingTypeZ):
    N_parameters = 2
    unstable     = False
    p_nlF_Hz     = 0
    p_nl_A       = 0
    max_amplitude = 1

    def update(self, nlF_Hz, nl_A):
        self.p_nlF_Hz = nlF_Hz
        self.p_nl_A = nl_A

    def reduce(self):
        return [
            self.p_nlF_Hz,
            self.p_nl_A,
        ]

    def option_set(self, minimum_BW_Hz = None, **kwargs):
        super(CodingnlFnlA, self).option_set(**kwargs)
        if minimum_BW_Hz is not None:
            self.max_amplitude = 1 - (minimum_BW_Hz / self.sys.F_nyquist_Hz)
        #TODO, should modify p_nl_A for the new max_amplitude
        return

    def update_roots(self, r1):
        """
        """
        #TODO should honor max_amplitude
        amp = abs(r1)
        F_Hz = np.angle(r1) / np.pi * self.sys.F_nyquist_Hz
        relF_Hz = F_Hz / self.sys.F_cutoff_Hz
        if relF_Hz >= 1:
            relF_Hz = .999
        elif relF_Hz <= -1:
            relF_Hz = .001
        self.p_nlF_Hz = relF_Hz / (1 - relF_Hz**2)**.5

        if amp < 1:
            #amp = x / (1. + x**2)**.5
            #amp**2 = x**2 / (1. + x**2)
            #amp**2 (1. + x**2) = x**2
            #amp**2  = x**2 ( 1 - amp**2)
            #amp**2 / (1 - amp**2) = x**2
            self.unstable = False
            self.p_nl_A = amp / (1 - amp**2)**.5
        else:
            #amp = (1. + x**2)**.5 / x
            #amp**2 =  (1. + x**2) / x**2
            #amp**2 x**2 = x**2 + 1
            #amp**2  = x**2 (amp**2 - 1)
            #amp**2 / (amp**2 - 1) = x**2
            self.unstable = True
            amp = 1 / amp
            self.p_nl_A = amp / (1 - amp**2)**.5
        return

    @property
    def gain_effect(self):
        #always 1 since it is complex
        if not self.unstable:
            return 1
        else:
            return 1

    @property
    def relF_Hz(self):
        val = self.p_nlF_Hz / (1. + self.p_nlF_Hz**2)**.5
        return val

    @property
    def amplitude(self):
        if not self.unstable:
            amp = self.max_amplitude * self.p_nl_A / (1. + self.p_nl_A**2)**.5
        else:
            amp = (1. + self.p_nl_A**2)**.5 / self.p_nl_A / self.max_amplitude
        return amp

    def transfer(self):
        #frequency, logarithmic amplitude
        amp  = self.amplitude
        r    = amp * np.cos(np.pi * self.sys.F_cutoff_Hz * self.relF_Hz / self.sys.F_nyquist_Hz)
        Xn   = self.sys.Xzn_grid
        Xnsq = self.sys.Xzn_grid_sq
        return (amp * amp) * Xnsq - 2 * Xn * r + 1

    def derivative(self):
        if self.disable:
            return []
        #real/imaginary part of root
        if not self.unstable:
            sqp5 = (1. + self.p_nl_A**2)**.5
            amp = self.max_amplitude * self.p_nl_A / sqp5
            DampDl = self.max_amplitude * (1 / sqp5 - self.p_nl_A**2 / (1. + self.p_nl_A**2)**1.5)
        else:
            sqp5 = (1. + self.p_nl_A**2)**.5
            amp = sqp5 / self.p_nl_A / self.max_amplitude
            DampDl = (1/sqp5 - sqp5 / self.p_nl_A**2) / self.max_amplitude

        Fsqp5     = (1. + self.p_nlF_Hz**2)**.5
        relF_Hz   = (self.p_nlF_Hz / Fsqp5)
        DrelFHzDp = (1 / Fsqp5 - self.p_nlF_Hz**2 / (1. + self.p_nlF_Hz**2)**1.5)

        F_nyquist_Hz = self.sys.F_nyquist_Hz
        Xn           = self.sys.Xzn_grid
        Xnsq         = self.sys.Xzn_grid_sq
        rcos     = np.cos(np.pi * self.sys.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
        rsin     = np.sin(np.pi * self.sys.F_cutoff_Hz * relF_Hz / F_nyquist_Hz)
        return [
            DrelFHzDp * (2 * amp * np.pi * self.sys.F_cutoff_Hz / self.sys.F_nyquist_Hz) * (Xn * rsin),
            ((2 * amp * DampDl) * Xnsq - Xn * (2 * rcos * DampDl)),
        ]

    def roots_c(self):
        #real/imaginary part of root
        return [self.amplitude * np.exp(Ipi * abs(self.relF_Hz) * self.sys.F_cutoff_Hz / self.sys.F_nyquist_Hz)]


