"""
"""

import numpy as np
import scipy
import scipy.signal
import scipy.linalg

from .eig_algorithms import eigspaces_right_real


def inverse_DSS(A, B, C, D, E):
    constrN = A.shape[0]
    statesN = A.shape[1]
    inputsN = B.shape[1]
    outputN = C.shape[0]

    constrNnew = constrN + inputsN
    statesNnew = statesN + inputsN
    inputsNnew = inputsN
    outputNnew = outputN

    assert(inputsN == outputN)

    newA = np.zeros((constrNnew, statesNnew))
    newE = np.zeros((constrNnew, statesNnew))
    newB = np.zeros((constrNnew, inputsNnew))
    newC = np.zeros((outputNnew, statesNnew))
    newD = np.zeros((outputNnew, inputsNnew))

    newA[:constrN, :statesN] = A
    newA[:constrN, statesN:] = B
    newA[constrN:, :statesN] = C
    newA[constrN:, statesN:] = D
    newE[:constrN, :statesN] = E
    newB[constrN:, :inputsN] = -1
    newC[:outputN, statesN:] = 1

    return newA, newB, newC, newD, newE


def reduce_modal(
        A, B, C, D, E,
        mode = 'C',
        tol = 1e-7,
        use_svd = False
):
    """
    TODO simplify a bunch!
    """
    v_pairs = eigspaces_right_real(A, E, tol = tol)
    A_projections = []
    if mode == 'C':
        #should maybe use SVD for rank estimation?

        #evects are columns are eig-idx, rows are eig-vectors
        for eigs, evects in v_pairs:
            #columns are B-in, rows are SS-eigen
            q, r, P = scipy.linalg.qr((E @ evects).T @ B, pivoting = True)
            for idx in range(q.shape[0] - 1, -1, -1):
                if np.sum(r[idx]**2) < tol:
                    continue
                else:
                    idx += 1
                    break
            idx_split = idx
            A_project = evects @ q[idx_split : ].T
            if A_project.shape[1] > 0:
                A_projections.append(A_project)
    elif mode == 'O':
        #TODO untested so far
        for eigs, evects in v_pairs:
            print(eigs)
            print(evects)
            #columns are C-out, rows are SS-eigen
            q, r, P = scipy.linalg.qr((C @ evects).T, pivoting = True)
            for idx in range(q.shape[0] - 1, -1, -1):
                if np.sum(r[idx]**2) < tol:
                    continue
                else:
                    idx += 1
                    break
            idx_split = idx
            A_project = evects @ q[idx_split : ].T
            if A_project.shape[1] > 0:
                A_projections.append(A_project)
    else:
        raise RuntimeError("Can only reduce mode='C' or 'O'")

    if not A_projections:
        return A, B, C, D, E, False

    A_projections = np.hstack(A_projections)
    Aq, Ar = scipy.linalg.qr(A_projections, mode = 'economic')
    A_SS_projection = np.diag(np.ones(A.shape[0])) - Aq @ Aq.T.conjugate()

    if use_svd:
        Au, As, Av = np.linalg.svd(A_SS_projection)
        for idx in range(len(As) - 1, -1, -1):
            if As[idx] < tol:
                continue
            else:
                idx += 1
                break
        idx_split = idx
        p_project_imU = Au[:, :idx_split]
        p_project_kerU = Au[:, idx_split :]
    else:
        Aq, Ar, Ap = scipy.linalg.qr(A_SS_projection, mode = 'economic', pivoting = True)
        for idx in range(Aq.shape[0] - 1, -1, -1):
            if np.sum(Ar[idx]**2) < tol:
                continue
            else:
                idx += 1
                break
        idx_split = idx
        p_project_imU = Aq[:, :idx_split]
        p_project_kerU = Aq[:, idx_split :]

    E_projections = E @ p_project_kerU
    Eq, Er = scipy.linalg.qr(E_projections, mode = 'economic')
    E_SS_projection = np.diag(np.ones(A.shape[0])) - Eq @ Eq.T.conjugate()

    if use_svd:
        Eu, Es, Ev = np.linalg.svd(E_SS_projection)
        p_project_im = Ev[:idx_split]
        p_project_ker = Ev[idx_split : ]
    else:
        Eq, Er, Ep = scipy.linalg.qr(E_SS_projection, mode = 'economic', pivoting = True)
        p_project_im = Eq.T[:idx_split]
        p_project_ker = Eq.T[idx_split : ]

    #TODO, have this check the mode argument
    if mode == 'C':
        assert(np.all((p_project_ker @ B)**2 < tol))
    if mode == 'O':
        assert(np.all((C @ p_project_kerU)**2 < tol))

    B = p_project_im @ B
    A = p_project_im @ A @ p_project_imU
    E = p_project_im @ E @ p_project_imU
    C = C @ p_project_imU
    return A, B, C, D, E, True


def controllable_staircase(
    A, B, C, D, E,
    tol = 1e-9,
):
    """
    Implementation of 
    COMPUTATION  OF IRREDUCIBLE  GENERALIZED STATE-SPACE REALIZATIONS ANDRAS VARGA
    using givens rotations.

    it is very slow, but numerically stable

    TODO, add pivoting,
    TODO, make it use the U-T property on E better for speed
    TODO, make it output Q and Z to apply to aux matrices, perhaps use them on C
    """
    #from icecream import ic
    #import tabulate
    Ninputs = B.shape[1]
    Nstates = A.shape[0]
    Nconstr = A.shape[1]
    Noutput = C.shape[0]

    BA, E = scipy.linalg.qr_multiply(
        E,
        np.hstack([B, A]),
        pivoting = False,
        mode = 'left'
    )

    Nmin = min(Nconstr, Nstates)
    for CidxBA in range(0, Nmin - 1):
        for RidxBA in range(Nconstr-1, CidxBA, -1):
            #create a givens rotation for Q reduction on BA
            BAv0 = BA[RidxBA - 1, CidxBA]
            BAv1 = BA[RidxBA, CidxBA]
            BAvSq = BAv0**2 + BAv1**2
            if BAvSq < tol:
                continue
            BAvAbs = BAvSq**0.5
            c = BAv1 / BAvAbs
            s = BAv0 / BAvAbs
            M = np.array([
                [s, +c],
                [-c, s]
            ])
            BA[RidxBA-1:RidxBA+1, :] = M @ BA[RidxBA-1:RidxBA+1, :]

            #TODO, use the U-T to be more efficient
            E[RidxBA-1:RidxBA+1, :] = M @ E[RidxBA-1:RidxBA+1, :]

            Cidx = RidxBA
            Ridx = RidxBA

            #row and col swap
            Ev0 = E[Ridx, Cidx-1]
            Ev1 = E[Ridx, Cidx]
            EvSq = Ev0**2 + Ev1**2
            if EvSq < tol:
                continue
            EvAbs = EvSq**0.5
            c = Ev0 / EvAbs
            s = Ev1 / EvAbs
            MT = np.array([
                [s, +c],
                [-c, s]
            ])
            BA[:, Ninputs:][:, Cidx-1:Cidx+1] = BA[:, Ninputs:][:, Cidx-1:Cidx+1] @ MT
            C[:, Cidx-1:Cidx+1] = C[:, Cidx-1:Cidx+1] @ MT
            #TODO, use the U-T to be more efficient
            E[:, Cidx-1:Cidx+1] = E[:, Cidx-1:Cidx+1] @ MT

    B = BA[:, :Ninputs]
    A = BA[:, Ninputs:]
    return A, B, C, D, E


