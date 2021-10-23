"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import copy
import declarative

from IIRrational.utilities.np import logspaced
from IIRrational.utilities.mpl import mplfigB
from IIRrational.statespace.dense import matrix_algorithms

import numpy.testing


def test_QRH():

    N = 10
    M = np.random.rand(N, N)
    M = np.array(
        [
            [1, 0, 0],
            [1, 1, 0],
            [0, 0, 0],
        ],
        float,
    )

    eye = np.eye(M.shape[0], M.shape[1])

    R, [Q], [QT] = matrix_algorithms.QR(
        mat=M,
        mshadow=None,
        qmul=[eye],
        qAmul=[eye],
        pivoting=False,
        # method   = 'Householder',
        method="Givens",
        Rexact=False,
    )
    R2, [Q], [QT] = matrix_algorithms.QR(
        mat=M,
        mshadow=None,
        qmul=[eye],
        qAmul=[eye],
        pivoting=False,
        # method   = 'Householder',
        method="Givens",
        Rexact=True,
    )

    import tabulate

    print("near", tabulate.tabulate(R))
    print("exact", tabulate.tabulate(R2))
    print(tabulate.tabulate(Q))
    print(tabulate.tabulate(QT))

    numpy.testing.assert_almost_equal(Q @ Q.T, eye)
    numpy.testing.assert_almost_equal(Q @ QT, eye)
    numpy.testing.assert_almost_equal(Q @ R2, M)


def test_QRHpivot():

    N = 10
    M = np.random.rand(N, N)
    M = np.array(
        [
            [1, 0, 0],
            [0, 1, 0],
            [0, 1, 0],
        ],
        float,
    )

    eye = np.eye(M.shape[0], M.shape[1])

    R, [Q], [QT], P = matrix_algorithms.QR(
        mat=M,
        mshadow=None,
        qmul=[eye],
        qAmul=[eye],
        pivoting=True,
        # method   = 'Householder',
        method="Givens",
        Rexact=True,
    )

    import tabulate

    print(P)
    print(tabulate.tabulate(R))
    # print(tabulate.tabulate(Q))
    # print(tabulate.tabulate(QT))

    numpy.testing.assert_almost_equal(Q @ Q.T, eye)
    numpy.testing.assert_almost_equal(Q @ QT, eye)
    numpy.testing.assert_almost_equal((Q @ R)[:, P], M)
