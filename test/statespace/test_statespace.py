"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import scipy
import scipy.signal

from IIRrational.utilities.np import logspaced
from IIRrational.utilities.mpl import mplfigB
from IIRrational.statespace import dense

import scipy.signal

c_m_s = 299792458


def print_ssd(ssd):
    print("B", ssd.B)
    print("A", ssd.A)
    print("E", ssd.E)
    print("C", ssd.C)
    print("D", ssd.D)


def settest(test_trigger, tpath_join, plot, Zc, Zr, Pc, Pr, k):
    Zc = 2 * np.pi * np.asarray(Zc)
    Zr = 2 * np.pi * np.asarray(Zr)
    Pc = 2 * np.pi * np.asarray(Pc)
    Pr = 2 * np.pi * np.asarray(Pr)

    Z = np.concatenate([Zc, Zc.conjugate(), Zr])
    P = np.concatenate([Pc, Pc.conjugate(), Pr])

    sys1 = dense.ZPKdict(
        name="zpk",
        zdict=dict(c=Zc, r=Zr),
        pdict=dict(c=Pc, r=Pr),
        k=k,
    )
    print_ssd(sys1)

    F_Hz = logspaced(0.01, 1e3, 1000)
    xfer = sys1.xfer(
        F_Hz=F_Hz,
        iname="zpk.i0",
        oname="zpk.o0",
    )

    # TODO, use IIRrational version since statespace is likely more numerically
    # stable than crappy scipy implementation
    w, xfer2 = scipy.signal.freqs_zpk(Z, P, k, 2 * np.pi * F_Hz)

    def trigger(fail, plot):
        axB = mplfigB(Nrows=2)
        axB.ax0.loglog(F_Hz, abs(xfer))
        axB.ax1.semilogx(F_Hz, np.angle(xfer, deg=True))
        axB.ax0.loglog(F_Hz, abs(xfer2))
        axB.ax1.semilogx(F_Hz, np.angle(xfer2, deg=True))
        axB.save(tpath_join("test"))

    with test_trigger(trigger, plot=plot):
        np.testing.assert_almost_equal(xfer, xfer2, decimal=5)

    sys_inv = sys1.inverse()
    xfer = 1 / sys_inv.xfer(
        F_Hz=F_Hz,
        iname="zpk.o0",
        oname="zpk.i0",
    )

    def trigger(fail, plot):
        axB = mplfigB(Nrows=2)
        axB.ax0.loglog(F_Hz, abs(xfer))
        axB.ax1.semilogx(F_Hz, np.angle(xfer, deg=True))
        axB.ax0.loglog(F_Hz, abs(xfer2))
        axB.ax1.semilogx(F_Hz, np.angle(xfer2, deg=True))
        axB.save(tpath_join("test_inv"))

    with test_trigger(trigger, plot=plot):
        np.testing.assert_almost_equal(xfer, xfer2, decimal=5)


def test_xfers(test_trigger, tpath_join, tpath_preclear, plot):
    settest(
        test_trigger,
        tpath_join,
        plot,
        Zc=[-1 + 1j, -1 + 5j],
        Zr=[-100, -200],
        Pc=[-1 + 2j, -1 + 6j],
        Pr=[-10, -20],
        # Zc = 2 * np.pi * np.asarray([])
        # Pc = 2 * np.pi * np.asarray([])
        k=1,
    )


def test_xfers2(test_trigger, tpath_join, tpath_preclear, plot):
    settest(
        test_trigger,
        tpath_join,
        plot,
        Zc=[-1 + 1j, -1 + 5j],
        Zr=[-100, -200],
        Pc=[-1 + 2j, -1 + 6j],
        Pr=[-10, -20, -30],
        k=1,
    )


def test_xfers12(test_trigger, tpath_join, tpath_preclear, plot):
    settest(test_trigger, tpath_join, plot, Zc=[], Zr=[-100], Pc=[-1 + 2j], Pr=[], k=1)


def test_xfers22(test_trigger, tpath_join, tpath_preclear, plot):
    settest(
        test_trigger, tpath_join, plot, Zc=[], Zr=[-10, -10], Pc=[-1 + 2j], Pr=[], k=1
    )


def test_xfers2i(test_trigger, tpath_join, tpath_preclear, plot):
    settest(
        test_trigger,
        tpath_join,
        plot,
        Zc=[-1 + 1j, -1 + 5j],
        Zr=[-100, -200, -10, -10],
        Pc=[-1 + 2j, -1 + 6j],
        Pr=[-10, -20, -30],
        k=1,
    )
