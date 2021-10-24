# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

# from wavestate.iirrational import fitters
from wavestate.iirrational.fitters_rational.polynomials import imchebychev
import numpy as np

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


def assert_roots_close(r1, r2, tol=1e-6):
    def key(v):
        return v.real, v.imag

    r1 = np.asarray(sorted(r1, key=key))
    r2 = np.asarray(sorted(r2, key=key))
    ratios = r1 / r2 - 1
    select = abs(ratios) > tol
    print(ratios)

    r1 = r1[select]
    r2 = r2[select]

    def key(v):
        return v.imag, v.real

    r1 = np.asarray(sorted(r1, key=key))
    r2 = np.asarray(sorted(r2, key=key))
    ratios = r1 / r2 - 1
    select = abs(ratios) > tol

    r1 = r1[select]
    r2 = r2[select]
    print(ratios)
    assert len(r1) == 0
    return


def test_imcheby_4roots_1scale():
    X_scale = 1
    roots = [-0.8 * X_scale, 0 + 0.5j * X_scale, 0 - 0.5j * X_scale, -0.5 * X_scale]
    c, lnG = imchebychev.fromroots_lnG(roots, X_scale=X_scale)

    F = np.linspace(0, X_scale, X_scale)
    X = 1j * F

    y_main, y_main_lnG = imchebychev.valfromroots_lnG(
        X, roots, X_scale=X_scale, lnG=lnG
    )
    y_poly, y_poly_lnG = imchebychev.val_lnG(X, c, X_scale=X_scale, lnG=lnG)
    V, VlnG = imchebychev.vander_lnG(X, len(c) - 1, X_scale=X_scale, lnG=lnG)
    y_vander = np.dot(V, c) * np.exp(VlnG)

    np.testing.assert_allclose(y_vander, y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_vander)
    return


# @pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_imcheby_4roots():
    X_scale = 100
    roots = [-0.8 * X_scale, 0 + 0.5j * X_scale, 0 - 0.5j * X_scale, -0.5 * X_scale]
    assert_roots_close(
        roots,
        imchebychev.roots(
            imchebychev.fromroots_lnG(roots, X_scale=X_scale)[0], X_scale=X_scale
        ),
    )

    c, lnG = imchebychev.fromroots_lnG(roots, X_scale=X_scale)
    print(c)

    F = np.linspace(0, X_scale, X_scale)
    X = 1j * F

    y_main, y_main_lnG = imchebychev.valfromroots_lnG(
        X, roots, X_scale=X_scale, lnG=lnG
    )
    y_poly, y_poly_lnG = imchebychev.val_lnG(X, c, X_scale=X_scale, lnG=lnG)
    V, VlnG = imchebychev.vander_lnG(X, len(c) - 1, X_scale=X_scale, lnG=lnG)
    y_vander = np.dot(V, c) * np.exp(VlnG)

    # from wavestate.iirrational.utilities.mpl import mplfigB
    # axB = mplfigB()
    # axB.ax0.plot(F, abs(y_main))
    # axB.ax0.plot(F, abs(y_poly))
    # axB.ax0.plot(F, abs(y_vander))

    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_vander, y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG) / y_vander, 1)
    return


def test_imcheby_3roots():
    X_scale = 100
    roots = [-0.8 * X_scale, 0 + 0.5j * X_scale, 0 - 0.5j * X_scale]
    assert_roots_close(
        roots,
        imchebychev.roots(
            imchebychev.fromroots_lnG(roots, X_scale=X_scale)[0], X_scale=X_scale
        ),
    )

    c, lnG = imchebychev.fromroots_lnG(roots, X_scale=X_scale)

    F = np.linspace(0, X_scale, X_scale)
    X = 1j * F

    y_main, y_main_lnG = imchebychev.valfromroots_lnG(
        X, roots, X_scale=X_scale, lnG=lnG
    )
    y_poly, y_poly_lnG = imchebychev.val_lnG(X, c, X_scale=X_scale, lnG=lnG)
    V, VlnG = imchebychev.vander_lnG(X, len(c) - 1, X_scale=X_scale, lnG=lnG)
    y_vander = np.dot(V, c) * np.exp(VlnG)

    # from wavestate.iirrational.utilities.mpl import mplfigB
    # axB = mplfigB()
    # axB.ax0.plot(F, abs(y_main))
    # axB.ax0.plot(F, abs(y_poly))
    # axB.ax0.plot(F, abs(y_vander))

    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_vander)
    return
