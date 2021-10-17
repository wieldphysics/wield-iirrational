"""
"""

import numpy as np
import scipy
import scipy.signal
from . import order_reduce

def ss2zpk(
    A, B, C, D, E = None,
    idx_in = None,
    idx_out = None,
    Q_rank_cutoff = 1e-5,
    Q_rank_cutoff_unstable = None,
    F_match_Hz = 1e-10,
    fmt = 'IIRrational',
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

    if E is None:
        p = scipy.linalg.eig(A, left = False, right = False)
    else:
        p = scipy.linalg.eig(A, E, left = False, right = False)
    SS = np.block([[A, B], [C, D]])
    if E is None:
        z = scipy.linalg.eig(
            a = SS,
            b = np.diag(np.concatenate([np.ones(A.shape[0]), np.zeros(1)])),
            left = False,
            right = False
        )
    else:
        SSE = np.block([
            [E,                                   np.zeros(E.shape[0]).reshape(-1, 1)],
            [np.zeros(E.shape[1]).reshape(1, -1), np.zeros(1).reshape(1, 1)]
        ])
        z = scipy.linalg.eig(
            a = SS,
            b = SSE,
            left = False,
            right = False
        )
    z = np.asarray([_ for _ in z if np.isfinite(_.real)])
    k = 1
    z, p, k = order_reduce.order_reduce_zpk(
        (z, p, k),
        reduce_c = True,
        reduce_r = True,
        Q_rank_cutoff = Q_rank_cutoff,
        Q_rank_cutoff_unstable = Q_rank_cutoff_unstable,
    )
    s_match_wHz = F_match_Hz * 2j * np.pi
    tf0 = (np.matmul(C, np.matmul(np.linalg.inv(np.eye(A.shape[0])*s_match_wHz - A), B)) + D)[..., 0, 0]
    w, zpk0 = scipy.signal.freqs_zpk(z, p, k, s_match_wHz)
    k = abs(tf0 / zpk0)

    if fmt == 'IIRrational':
        z = np.asarray(z) / (2 * np.pi)
        p = np.asarray(p) / (2 * np.pi)
        k = np.asarray(k)  * (2 * np.pi)**(len(z) - len(p))
    elif fmt == 'scipy':
        pass
    else:
        raise RuntimeError("Unrecognized fmt parameter")
    return z, p, k


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
