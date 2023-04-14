# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import numpy as np
from wield.iirrational.representations.polynomials import chebychev

from numpy.polynomial import chebyshev

chebroots = chebyshev.chebroots
chebfromroots = chebyshev.chebfromroots
chebder = chebyshev.chebder
chebadd = chebyshev.chebadd
chebmul = chebyshev.chebmul
chebdiv = chebyshev.chebdiv
chebval = chebyshev.chebval


def abs_sq(v):
    return np.real(v) ** 2 + np.imag(v) ** 2


test_roots = [
    0.5
    * np.asarray(
        (
            -2,
            -1.1 + 1j,
            -1.1 + 1j,
            -1.1 + 1j,
            -2 + 1j,
            -1 + 1j,
            -1 + 1j,
            -1 + 1j,
            -1.1 - 1j,
            -1.1 - 1j,
            -1.1 - 1j,
            -2 + 1j,
            -1 - 1j,
            -1 - 1j,
            -1 - 1j,
        )
    ),
]


@pytest.mark.parametrize("roots", test_roots)
def test_root_chaser(roots):
    c = chebfromroots(roots)
    c_orig = c

    def calc(z, pr=False):
        cval = chebval(z, c)
        cdval = chebval(z, cd)
        rat1 = cval / cdval
        ret = None
        if np.any(abs_sq(cdval) > tol):
            ret = rat1
        cd2val = chebval(z, cd2)
        rat2 = cdval / cd2val
        if ret is None and np.any(abs_sq(cd2val) > tol):
            ret = rat2
        cd3val = chebval(z, cd3)
        rat3 = cd2val / cd3val
        if ret is None and np.any(abs_sq(cd3val) > tol):
            ret = rat3
        cd4val = chebval(z, cd4)
        rat4 = cd3val / cd4val
        if ret is None:
            ret = rat4
        ret = rat2 * 1
        if pr:
            # print(abs(cval), abs(cdval), abs(cd2val), abs(cd3val))
            print(abs(rat1), abs(rat2), abs(rat3), abs(rat4))
        return ret

    def laguerre(z):
        cval = chebval(z, c)
        cdval = chebval(z, cd)
        cd2val = chebval(z, cd2)

        G = cdval / cval
        Gsq = G ** 2
        H = Gsq - cd2val / cval
        k = 2
        m = 2
        n = m + k

        disc = k / m * (n * H - Gsq)
        disc_root = disc ** 0.5
        a_p = n / (G + disc_root)
        a_n = n / (G - disc_root)
        if abs(a_p) < abs(a_n):
            return a_p
        else:
            return a_n

    print()

    tol = 1e-8
    current_z = (-1.01 + 1.01j) / 2
    roots = chebroots(c)
    select = np.nonzero(abs(roots - current_z) < 0.05)[0]
    print(select[1:])
    idx_update = 0

    for n in range(1):
        for idx_update in range(len(select)):
            myroot_idx = select[idx_update]
            subroots = np.concatenate(
                [roots[select[:idx_update]], roots[select[idx_update + 1 :]]]
            )
            print(subroots)
            # subroots = [
            #    -1 + 1j,
            #    -1 + 1j,
            # ]
            c = c_orig
            # for sr in subroots:
            #    sc = chebfromroots([sr])
            #    q, r = chebdiv(c, sc)
            #    print('remainder', abs(r), r)
            #    c = q
            # rescale by the root density nearest to the origin
            c = c / c[1]
            cd = chebder(c)
            cd2 = chebder(cd)
            cd3 = chebder(cd2)
            cd4 = chebder(cd3)

            current_z = (-1.01 + 1.01j) / 2
            current_z = -1.1
            scoots = []
            for idx in range(100):
                rat = calc(current_z, pr=True)
                # print("CVAL: ", abs(rat))
                # diff_z = -scale * rat
                # new_z = current_z + diff_z
                # val = calc(new_z)
                # a, b = np.polyfit(diff_z, val, deg = 1)
                # scoot = -b/a
                # scoot = scoot * abs(rat) / (abs(scoot)/100 + abs(rat))
                # scoot = -rat
                scoot = -laguerre(current_z)
                # scoot = -rat
                scoots.append(scoot)
                # print('scoot', scoot)
                current_z = current_z + scoot
            scoots = np.asarray(scoots)
            print(abs(scoots))
            roots[myroot_idx] = current_z
            print(myroot_idx, current_z)
            # print(abs(scoots))

    # assert(False)
