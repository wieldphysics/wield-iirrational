"""
Utilities to manipulate ZPK roots S to/from Z, and make transfer functions
"""

import numpy as np
import declarative


def analytic_translation(
        F_Hz,
        time_s = None,
        F_displacement_Hz = 0,
):
    """
    time_bins shou
    """
    if time_s is None:
        time_s = np.linspace(0, (len(F_Hz) - 1)/np.max(F_Hz) / 2, 3 * (1 + 1 * len(F_Hz) / 2)) * 1.2
    if isinstance(time_s, (float, int)):
        time_s = np.linspace(0, (len(F_Hz) - 1)/np.max(F_Hz) / 2, time_s * (1 + 1 * len(F_Hz) / 2)) * 1.2

    max_F_Hz = np.max(F_Hz)

    t, f = np.meshgrid(time_s, F_Hz)
    fourier_c = np.cos(2 * np.pi * f * t)
    fourier_s = np.sin(2 * np.pi * f * t)

    disp_factor = np.exp(-F_displacement_Hz * time_s).reshape(-1, 1)
    #u_c, s_c, v_c = np.linalg.svd(fourier_c)
    #u_s, s_s, v_s = np.linalg.svd(fourier_s)
    fourier_c_inv = disp_factor * np.linalg.pinv(fourier_c)
    fourier_s_inv = disp_factor * np.linalg.pinv(fourier_s)

    trans_c = np.dot(fourier_c, fourier_c_inv)
    trans_s = np.dot(fourier_s, fourier_s_inv)

    trans_sc = np.dot(fourier_s, fourier_c_inv)
    trans_cs = np.dot(fourier_c, fourier_s_inv)

    def apply_translation(cplx_vect):
        return (np.dot(trans_c, cplx_vect.real) + 1j * np.dot(trans_s, cplx_vect.imag))

    def apply_continuation(cplx_vect):
        return -(1j * np.dot(trans_sc, cplx_vect.real) + np.dot(trans_cs, cplx_vect.imag))

    def stability_test(cplx_vect):
        return apply_translation(cplx_vect) / apply_continuation(cplx_vect)

    return wavestate.bunch.Bunch(locals())


def analytic_translation_dual(
        F_Hz,
        time_s = None,
        F_displacement_Hz = 0,
):
    """
    time_bins shou
    """
    if time_s is None:
        time_s = np.linspace(0, (len(F_Hz) - 1)/np.max(F_Hz) / 2, 3 * (1 + 1 * len(F_Hz) / 2)) * 1.2
    if isinstance(time_s, (float, int)):
        time_s = np.linspace(0, (len(F_Hz) - 1)/np.max(F_Hz) / 2, time_s * (1 + 1 * len(F_Hz) / 2)) * 1.2

    max_F_Hz = np.max(F_Hz)

    t, f = np.meshgrid(time_s, F_Hz)
    fourier_c = np.cos(2 * np.pi * f * t)
    fourier_s = np.sin(2 * np.pi * f * t)
    fourier_cs = np.hstack([fourier_c, fourier_s])
    fourier_sc = np.hstack([fourier_s, -fourier_c])
    print(fourier_cs.shape)

    disp_factor = np.exp(-F_displacement_Hz * time_s).reshape(-1, 1)
    print(disp_factor.shape)
    disp_factor = np.vstack([disp_factor, disp_factor])
    print(disp_factor.shape)

    fourier_cs_inv = disp_factor * np.linalg.pinv(fourier_cs)

    trans    = np.dot(fourier_cs, fourier_cs_inv)
    trans_sc = np.dot(fourier_sc, fourier_cs_inv)

    def apply_translation(cplx_vect):
        return np.dot(trans, cplx_vect)

    def apply_continuation(cplx_vect):
        return -1j * np.dot(trans_sc, cplx_vect)

    return wavestate.bunch.Bunch(locals())
