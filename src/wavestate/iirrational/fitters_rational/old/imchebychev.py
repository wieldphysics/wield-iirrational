# -*- coding: utf-8 -*-
"""
These are chebychev polynomials on a rotated domain. While they should be the chebychev polynomials form
F(ix) = Sum_i c_i * C_i(x / i)

so that evaluated along the imaginary line they are standard chebychev's, these are not those. Instead they are the chebychevs of

F(ix) = Sum_i c_i * (C_ei(x) + iC_oi(x))

such that all of the c_i coefficients are pure real, but the even polynomials provide the real part and the odd ones the imaginary part. This is easy to create a vandermonde matrix from and easy to solve. The roots of this form of polynomial will be paired left/right over the imaginary line rather than paired with conjugates as for pure real coefficients.

To relate the coefficients to the original requires them to be rotated 90 degrees
"""


import numpy as np
from . import chebychev
#from . import standard

valfromroots_lnG            = chebychev.valfromroots_lnG
coeff_canonicalization_gain = chebychev.coeff_canonicalization_gain

def roots(c, X_scale = 1):
    c = np.asarray(c, complex)
    c2 = np.copy(c)
    if len(c) % 2 == 1:
        c2[1::2] *= 1j
    else:
        c2[0::2] *= -1j
    roots = chebychev.roots(c2, X_scale = X_scale)

    SGN = (c[0] / c[-1])
    count = roots.imag < 0 
    print(c)
    print("CHECK COUNT!", (c[0]), SGN, np.sum(count), len(roots))

    #for root in roots:
    #    print(root, abs(chebychev.val_lnG(root, c2, X_scale = X_scale)[0]))

    roots = 1j * roots

    return roots

def roots_lnG(c, X_scale = 1):
    c = np.asarray(c, complex)
    c2 = np.copy(c)
    if len(c) % 2 == 1:
        c2[1::2] *= 1j
    else:
        c2[0::2] *= -1j
    roots, lnG = chebychev.roots_lnG(c2, X_scale = X_scale)
    roots = 1j * roots

    #TODO, repair the real poly conjugate matching
    return roots, lnG

def fromroots_lnG(roots, X_scale = 1):
    roots = np.asarray(roots, complex)
    roots2 = roots * -1j
    c, lnG = chebychev.fromroots_lnG(roots2, X_scale = X_scale)
    if len(roots) % 2 == 0:
        c[1::2] *= -1j
    else:
        c[0::2] *= +1j
    assert(np.all(c.imag == 0))
    return c.real, lnG

def val_lnG(X, c, X_scale = 1, lnG = 0):
    c2 = np.asarray(c, complex)
    c2 = np.copy(c2)
    c2[1::2] *= +1j
    if len(c2) % 2 == 1:
        #print('A2')
        return chebychev.val_lnG(X * -1j, c2, X_scale = X_scale, lnG = lnG)
    else:
        #print('B2')
        val, lnG = chebychev.val_lnG(X * -1j, c2, X_scale = X_scale, lnG = lnG)
        return -val, lnG

def vander_lnG(X, N, X_scale = 1, lnG = 0):
    V, lnG = chebychev.vander_lnG(X * +1j, N, X_scale = X_scale, lnG = lnG)
    V = V.astype(complex)
    if N % 2 == 1:
        V[:, 1::2] *= 1j
        V[:, 0::2] *= -1
    else:
        V[:, 1::2] *= -1j
    return V, lnG

def companion(c):
    c2 = np.asarray(c, complex)
    c2 = np.copy(c2)
    if len(c) % 2 == 1:
        c2[1::2] *= 1j
    else:
        c2[0::2] *= 1j
    return 1j * np.polynomial.chebyshev.chebcompanion(c2)
