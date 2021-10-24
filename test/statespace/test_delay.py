"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import scipy
import scipy.signal

from wavestate.iirrational.utilities.np import logspaced
from wavestate.iirrational.utilities.mpl import mplfigB
from wavestate.iirrational.statespace import dense, StateSpaceDense

import scipy.signal

c_m_s = 299792458


def print_ssd(ssd):
    print("B", ssd.B)
    print("A", ssd.A)
    print("E", ssd.E)
    print("C", ssd.C)
    print("D", ssd.D)


def test_delay(tpath_join, test_trigger):
    length_m = 3995
    delta_t = length_m / c_m_s
    delta_t = 1
    axB = mplfigB(Nrows=2)

    for idx_ord in range(1, 7):
        arm1 = dense.delay("arm1", delta_t, order=idx_ord, method="bessel")
        print_ssd(arm1)

        F_Hz = logspaced(0.01 / delta_t, 2 / delta_t, 1000)
        xfer = arm1.xfer(
            F_Hz=F_Hz,
            iname="arm1.i0",
            oname="arm1.o0",
        )

        axB.ax0.semilogx(F_Hz, abs(xfer), label="order {}".format(idx_ord))
        axB.ax1.plot(F_Hz, np.angle(xfer, deg=True))

    xfer_delay = np.exp(-2j * np.pi * F_Hz * delta_t)
    axB.ax1.plot(F_Hz, np.angle(xfer_delay, deg=True), color="magenta", ls="--")
    axB.ax1.axvline(1 / delta_t / 4)
    axB.ax1.axvline(2 / delta_t / 4)
    axB.ax1.axvline(3 / delta_t / 4)
    axB.ax1.axvline(4 / delta_t / 4)
    axB.ax0.legend()
    axB.save(tpath_join("test"))


def test_big_delay(tpath_join, test_trigger):
    length_m = 3995
    delta_t = length_m / c_m_s
    delta_t = 1
    axB = mplfigB(Nrows=2)

    idx_ord = 100
    arm1 = dense.delay("arm1", delta_t, order=idx_ord, method="bessel")
    print_ssd(arm1)

    F_Hz = np.linspace(0.00 / delta_t, 50 / delta_t, 1000)
    xfer = arm1.xfer(
        F_Hz=F_Hz,
        iname="arm1.i0",
        oname="arm1.o0",
    )

    axB.ax0.semilogx(F_Hz, abs(xfer), label="order {}".format(idx_ord))
    axB.ax1.plot(F_Hz, np.angle(xfer, deg=True))

    xfer_delay = np.exp(-2j * np.pi * F_Hz * delta_t)
    axB.ax1.plot(F_Hz, np.angle(xfer_delay, deg=True), color="magenta", ls="--")
    axB.ax1.axvline(1 / delta_t / 4)
    axB.ax1.axvline(2 / delta_t / 4)
    axB.ax1.axvline(3 / delta_t / 4)
    axB.ax1.axvline(4 / delta_t / 4)
    axB.ax0.legend()
    axB.save(tpath_join("test"))
