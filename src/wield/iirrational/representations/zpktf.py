#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Calculation of the RMS of a filter with unit white noise passing through a ZPK.
This is done using an integral and residue calculus. The filter must have more
P1 than Z, unless the ZPK is in the Z domain,
where there is a natural cutoff frequency.
"""

import numpy as np

import numbers

from wield.bunch.depbunch import (
    DepBunch,
    depB_property,
    # NOARG,
)

from .root_bunch import (
    RBAlgorithms,
    RootBunch,
    root_constraints,
)

from .ratmath import (
    ZPKsum,
    ZPKprod,
    ZPKscalarprod,
    ZPKdiv,
    ZPKscalardiv,
    ZPKscalarsum,
    ZPKdivscalar,
)


class ZPKTF(DepBunch):
    RBalgo = RBAlgorithms()
    root_constraint = RBalgo.root_constraints.mirror_real

    def __build__(
        self,
        ZPK=None,
        zeros=None,
        poles=None,
        gain=None,
        F_nyquist_Hz="unknown",
        **kwargs
    ):
        super(ZPKTF, self).__build__(**kwargs)
        if zeros is None and ZPK is not None:
            try:
                zeros = ZPK.zeros
            except AttributeError:
                zeros = ZPK[0]
        if zeros is None:
            zeros = ()
        if poles is None and ZPK is not None:
            try:
                poles = ZPK.poles
            except AttributeError:
                poles = ZPK[1]
        if poles is None:
            poles = ()
        if gain is None and ZPK is not None:
            try:
                gain = ZPK.gain
            except AttributeError:
                gain = ZPK[2]
        if gain is None:
            gain = 1
        self.zeros = zeros
        self.poles = poles
        self.gain = gain
        self.F_nyquist_Hz = F_nyquist_Hz

    def __iter__(self):
        yield self.zeros.fullplane
        yield self.poles.fullplane
        yield self.gain

    def __str__(self):
        return "ZPKTF(Z={},P={},K={})".format(
            str(self.zeros), str(self.poles), str(self.gain)
        )

    @depB_property
    def test(self, val=1):
        return val

    @depB_property
    def gain(self, val):
        return val

    @depB_property
    def zeros(self, val):
        val = self.RBalgo.expect_atleast(val, constraint=self.root_constraint)
        return val

    @depB_property
    def poles(self, val):
        val = self.RBalgo.expect_atleast(val, constraint=self.root_constraint)
        return val

    @depB_property
    def order(self):
        return max(
            len(self.poles),
            len(self.zeros),
        )

    @depB_property
    def order_sos(self):
        return (
            max(
                len(self.poles),
                len(self.zeros),
            )
            + 1
        ) // 2

    @depB_property
    def order_relative(self):
        return len(self.zeros) - len(self.poles)

    @depB_property
    def order_total(self):
        return len(self.poles) + len(self.zeros)

    def xfer_eval(self, F_Hz):
        # TODO must add pole-zero rephasing for Z filters

        # TODO, make this name consistent in all classes
        if self.F_nyquist_Hz is None:
            X_grid = 1j * F_Hz
        else:
            # use Z^-1
            X_grid = np.exp(1j * np.pi * F_Hz / self.F_nyquist_Hz)
        h, lnG = self.poles.val_lnG(X_grid)
        h, lnG = self.zeros.val_lnG(X_grid, h=1 / h, lnG=-lnG)
        return h * (np.exp(lnG) * self.gain)

    def __add__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalarsum(self, other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, tuple):
            F_nyquist_Hz = common_nyquist(self, other)
            Z3, P3, K3 = ZPKsum(self, other)
        else:
            return NotImplemented
        return self.__class__(Z3, P3, K3, F_nyquist_Hz=F_nyquist_Hz)

    def __radd__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalarsum(self, other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, tuple):
            F_nyquist_Hz = common_nyquist(other, self)
            Z3, P3, K3 = ZPKsum(self, other)
        else:
            return NotImplemented
        return self.__class__(Z3, P3, K3, F_nyquist_Hz=F_nyquist_Hz)

    def __sub__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalarsum(self, -other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, tuple):
            F_nyquist_Hz = common_nyquist(self, other)
            Z2, P2, K2 = other
            Z3, P3, K3 = ZPKsum(self, (Z2, P2, -K2))
        else:
            return NotImplemented
        return self.__class__(Z3, P3, K3, F_nyquist_Hz=F_nyquist_Hz)

    def __rsub__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalarsum(self, -other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, tuple):
            F_nyquist_Hz = common_nyquist(other, self)
            Z2, P2, K2 = other
            Z3, P3, K3 = ZPKsum(self, (Z2, P2, -K2))
        else:
            return NotImplemented
        return self.__class__(Z3, P3, -K3, F_nyquist_Hz=F_nyquist_Hz)

    def __mul__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalarprod(self, other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, (tuple, ZPKTF)):
            F_nyquist_Hz = common_nyquist(self, other)
            Z3, P3, K3 = ZPKprod(self, other)
        else:
            return NotImplemented
        return self.__class__(zeros=Z3, poles=P3, gain=K3, F_nyquist_Hz=F_nyquist_Hz)

    def __rmul__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalarprod(self, other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, (tuple, ZPKTF)):
            F_nyquist_Hz = common_nyquist(other, self)
            Z3, P3, K3 = ZPKprod(other, self)
        else:
            return NotImplemented
        return self.__class__(zeros=Z3, poles=P3, gain=K3, F_nyquist_Hz=F_nyquist_Hz)

    def __truediv__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKscalardiv(self, other)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, (tuple, ZPKTF)):
            F_nyquist_Hz = common_nyquist(self, other)
            Z3, P3, K3 = ZPKdiv(self, other)
        else:
            return NotImplemented
        return self.__class__(zeros=Z3, poles=P3, gain=K3, F_nyquist_Hz=F_nyquist_Hz)

    def __rtruediv__(self, other):
        if isinstance(other, numbers.Real):
            Z3, P3, K3 = ZPKdivscalar(other, self)
            F_nyquist_Hz = self.F_nyquist_Hz
        elif isinstance(other, (tuple, ZPKTF)):
            F_nyquist_Hz = common_nyquist(other, self)
            Z3, P3, K3 = ZPKdiv(other, self)
        else:
            return NotImplemented
        return self.__class__(zeros=Z3, poles=P3, gain=K3, F_nyquist_Hz=F_nyquist_Hz)

    # TODO, version check here?
    __div__ = __truediv__
    __rdiv__ = __rtruediv__

    def __pow__(self, other):
        if isinstance(other, numbers.Complex):
            if other.imag != 0:
                return NotImplemented
            other = other.real

        if not isinstance(other, numbers.Integral):
            residual = other % 1
            if residual < -1e-14 or residual > 1e-14:
                return NotImplemented
            other = int((other + 1e-14) // 1)

        if other == 0:
            return self.__class__((), (), 1)
        elif other == 1:
            return self
        elif other == -1:
            return self.__class__(
                self.poles, self.zeros, 1 / self.gain, F_nyquist_Hz=self.F_nyquist_Hz
            )
        elif other > 1:
            return self.__class__(
                tuple(self.zeros) * other,
                tuple(self.poles) * other,
                self.gain ** other,
                F_nyquist_Hz=self.F_nyquist_Hz,
            )
        elif other < 1:
            other = -other
            return self.__class__(
                tuple(self.poles) * other,
                tuple(self.zeros) * other,
                self.gain ** (-other),
                F_nyquist_Hz=self.F_nyquist_Hz,
            )

    def __neg__(self):
        return self.__class__(
            self.poles, self.zeros, -self.gain, F_nyquist_Hz=self.F_nyquist_Hz
        )

    def __pos__(self):
        return self

    def abs_sq(self, F_nyquist_Hz="unknown"):
        if self.F_nyquist_Hz != "unknown":
            F_nyquist_Hz = self.F_nyquist_Hz
        if F_nyquist_Hz == "unknown":
            raise RuntimeError(
                " F_nyquist_Hz must be known to perform the correct"
                " ZPK manipulation for abs_sq"
            )
        elif F_nyquist_Hz is None:
            return self.__class__(
                tuple(self.poles) + tuple(np.asarray(self.poles).conjugate()),
                tuple(self.zeros) + tuple(np.asarray(self.zeros).conjugate()),
                self.gain ** 2,
                F_nyquist_Hz=self.F_nyquist_Hz,
            )
        else:
            return self.__class__(
                tuple(self.poles) + tuple(1 / np.asarray(self.poles)),
                tuple(self.zeros) + tuple(1 / np.asarray(self.zeros)),
                self.gain ** 2,
                F_nyquist_Hz=self.F_nyquist_Hz,
            )

    def assert_F_nyquist_Hz(self, F_nyquist_Hz):
        if self.F_nyquist_Hz == "unknown":
            return
        if self.F_nyquist_Hz != F_nyquist_Hz:
            raise RuntimeError("Incompatible ZPK representations!")
        return F_nyquist_Hz


def common_nyquist(LHS, RHS):
    """
    Assumes LHS and RHS are a ZPK or has a
    F_nyquist_Hz property.
    """
    if isinstance(LHS, ZPKTF):
        LHS_F_nyquist_Hz = LHS.F_nyquist_Hz
    else:
        LHS_F_nyquist_Hz = "unknown"

    if isinstance(RHS, ZPKTF):
        RHS_F_nyquist_Hz = RHS.F_nyquist_Hz
    else:
        RHS_F_nyquist_Hz = "unknown"

    if LHS_F_nyquist_Hz == "unknown":
        return RHS_F_nyquist_Hz
    elif RHS_F_nyquist_Hz == "unknown":
        return LHS_F_nyquist_Hz
    if LHS_F_nyquist_Hz != RHS_F_nyquist_Hz:
        raise RuntimeError(
            (
                "Math Operation should not be done on ZPKs arising from "
                "different representations, LHS is {}, RHS is {}"
            ).format(LHS_F_nyquist_Hz, RHS_F_nyquist_Hz)
        )
    return LHS_F_nyquist_Hz


# TODO,
def asZPKTF(
    ZPK,
    complete=False,
    F_nyquist_Hz=None,
    delay_s=None,
):
    from .zpk_with_data import ZPKwData

    if isinstance(ZPK, ZPKTF):
        return ZPK
    elif isinstance(ZPK, (tuple, list)):
        Z, P, K = ZPK
        return ZPKTF(
            zeros=asMRRB(Z, complete=complete),
            poles=asMRRB(P, complete=complete),
            gain=K,
            F_nyquist_Hz=F_nyquist_Hz,
        )
    elif isinstance(ZPK, ZPKwData):
        return ZPKTF(
            zeros=ZPK.zeros * ZPK.zeros_overlay,
            poles=ZPK.poles * ZPK.poles_overlay,
            gain=ZPK.gain,
            F_nyquist_Hz=F_nyquist_Hz,
        )
    # last ditch effort if it is some other wield.iirrational type
    return asZPKTF(ZPK.ZPKrep, complete=False)


def asMRRB(
    roots=None,
    r=(),
    c=(),
    complete=False,
):
    """
    Convenience Method to generate root_bunches with mirror real constraints
    from raw root lists
    """
    rb0 = RootBunch(r=r, c=c, constraint=root_constraints.mirror_real)

    if roots is not None:
        roots = np.asarray(roots)
        select_real = roots.imag == 0
        rr = roots[select_real]
        rc = roots[~select_real]
        if np.all(rc.imag < 0) and complete:
            return RootBunch(
                r=np.concatenate([rb0.r, rr]),
                c=np.concatenate([rb0.c, rc.conjugate()]),
                constraint=root_constraints.mirror_real,
            )
        elif np.all(rc.imag > 0) and complete:
            return RootBunch(
                r=np.concatenate([rb0.r, rr]),
                c=np.concatenate([rb0.c, rc]),
                constraint=root_constraints.mirror_real,
            )
        else:
            rb = ZPKTF.RBalgo.expect(roots, constraint=root_constraints.mirror_real)
            return rb * rb0
    else:
        return rb0
