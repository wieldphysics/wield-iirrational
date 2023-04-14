# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

# from wield.iirrational import fitters
from wield.iirrational.fitters_rational.polynomials import standard
import numpy as np

try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


def test_std_4roots_1scale():
    X_scale = 1
    roots = [-0.8 * X_scale, 0 + 0.5j * X_scale, 0 - 0.5j * X_scale, -0.5 * X_scale]
    c, lnG = standard.fromroots_lnG(roots, X_scale=X_scale)

    F = np.linspace(0, X_scale, X_scale)
    X = 1j * F

    y_main, y_main_lnG = standard.valfromroots_lnG(X, roots, X_scale=X_scale)
    y_poly, y_poly_lnG = standard.val_lnG(X, c, X_scale=X_scale, lnG=lnG)
    V, VlnG = standard.vander_lnG(X, len(c) - 1, X_scale=X_scale, lnG=lnG)
    y_vander = np.dot(V, c) * np.exp(VlnG)

    np.testing.assert_allclose(y_vander, y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_vander)
    return


# @pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
def test_std_4roots():
    X_scale = 100
    roots = [-0.8 * X_scale, 0 + 0.5j * X_scale, 0 - 0.5j * X_scale, -0.5 * X_scale]
    # print(standard.roots(standard.fromroots_lnG(roots, X_scale = X_scale)[0], X_scale = X_scale))

    c, lnG = standard.fromroots_lnG(roots, X_scale=X_scale)
    print(c, np.exp(lnG))

    F = np.linspace(0, X_scale, X_scale)
    X = 1j * F

    y_main, y_main_lnG = standard.valfromroots_lnG(X, roots, X_scale=X_scale, lnG=lnG)
    y_poly, y_poly_lnG = standard.val_lnG(X, c, X_scale=X_scale, lnG=lnG)
    V, VlnG = standard.vander_lnG(X, len(c) - 1, X_scale=X_scale, lnG=lnG)
    y_vander = np.dot(V, c) * np.exp(VlnG)

    # from wield.iirrational.utilities.mpl import mplfigB
    # axB = mplfigB()
    # axB.ax0.plot(F, abs(y_main))
    # axB.ax0.plot(F, abs(y_poly))
    # axB.ax0.plot(F, abs(y_vander))

    np.testing.assert_allclose(y_vander, y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG) / y_vander, 1)
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_poly * np.exp(y_poly_lnG))
    return


def test_std_3roots():
    X_scale = 100
    roots = [-0.8 * X_scale, 0 + 0.5j * X_scale, 0 - 0.5j * X_scale]
    print(
        standard.roots(
            standard.fromroots_lnG(roots, X_scale=X_scale)[0], X_scale=X_scale
        )
    )

    c, lnG = standard.fromroots_lnG(roots, X_scale=X_scale)

    F = np.linspace(0, X_scale, X_scale)
    X = 1j * F

    y_main, y_main_lnG = standard.valfromroots_lnG(X, roots, X_scale=X_scale, lnG=lnG)
    y_poly, y_poly_lnG = standard.val_lnG(X, c, X_scale=X_scale, lnG=lnG)
    V, VlnG = standard.vander_lnG(X, len(c) - 1, X_scale=X_scale, lnG=lnG)
    y_vander = np.dot(V, c) * np.exp(VlnG)

    # from wield.iirrational.utilities.mpl import mplfigB
    # axB = mplfigB()
    # axB.ax0.plot(F, abs(y_main))
    # axB.ax0.plot(F, abs(y_poly))
    # axB.ax0.plot(F, abs(y_vander))

    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_poly * np.exp(y_poly_lnG))
    np.testing.assert_allclose(y_main * np.exp(y_main_lnG), y_vander)
    return
