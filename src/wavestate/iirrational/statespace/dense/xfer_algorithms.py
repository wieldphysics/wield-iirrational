"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np


def ss2xfer(
    A, B, C, D, E = None,
    F_Hz          = None,
    idx_in        = None,
    idx_out       = None
):
    if idx_in is None:
        if B.shape[1] == 1:
            idx_in = 0
        else:
            raise RuntimeError("Must specify idx_in if B indicates MISO/MIMO system")
    if idx_out is None:
        if C.shape[0] == 1:
            idx_out = 0
        else:
            raise RuntimeError("Must specify idx_in if C indicates SIMO/MIMO system")
    B = B[:, idx_in:idx_in+1]
    C = C[idx_out:idx_out+1, :]
    D = D[idx_out:idx_out+1, idx_in:idx_in+1]

    s = 2j * np.pi * F_Hz
    if E is None:
        S = np.eye(A.shape[0]).reshape(1, *A.shape) * s.reshape(-1, 1, 1)
        return (np.matmul(C, np.matmul(np.linalg.inv(S - A.reshape(1, *A.shape)), B)) + D)[..., 0, 0]
    else:
        return (
            np.matmul(
                C,
                np.matmul(
                    np.linalg.inv(E*s.reshape(-1, 1, 1) - A.reshape(1, *A.shape)),
                    B
                )
            ) + D
        )[..., 0, 0]
