"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import scipy
import scipy.signal

from IIRrational.utilities.np import logspaced
from IIRrational.utilities.mpl import mplfigB
from IIRrational.TFmath import statespace
import IIRrational.statespace.dense as SS
from os import path

#def test_2x2(tpath):
#    M_zpk = [
#        [
#            #((),      (-3, -2), 200),
#            ((-1,    -5),      (-1+5j, -1-5j, -3), 200),
#            ((-1+1j, -1-1j, -5), (-1+5j, -1-5j,  -1+1j, -1-1j, -3),  1),
#        ], [
#            ((-1+1j, -1-1j, -5), (-1+5j, -1-5j,  -3),  1),
#            ((-5,),      (-1+1j, -1-1j, -3), 1),
#        ],
#    ]

def test_ABCDE(tpath):
    print()
    F_Hz = logspaced(.1, 100, 1000)

    z, p, k = ((-1,    -5),      (-1+5j, -1-5j, -3), 200)
    z = np.asarray(z)
    p = np.asarray(p)

    A, B, C, D, E = SS.zpk2cDSS(z, p, k, mode = 'CCF', rescale = 3j)
    A, B, C, D, E = SS.DSS_c2r(A, B, C, D, E)
    print("A", A)
    print("B", B)
    print("C", C)
    print("D", D)
    print("E", E)
    reduced = True
    while reduced:
        print("REDUCE")
        A, B, C, D, E, reduced = SS.reduce_modal(A, B, C, D, E, mode = 'O')
        #if not reduced:
        #    break
        #A, B, C, D, E, reduced = reduce_modal(A, B, C, D, E, mode = 'C')

    w, vr = scipy.linalg.eig(A, E, left = False, right = True)
    print("EIGS", w)
    print("A", A)
    print("E", E)
    print("B", B)

    axB = mplfigB(Nrows = 2)
    w, h = scipy.signal.freqs_zpk(
        z * 2 * np.pi,
        p * 2 * np.pi,
        k,
        F_Hz * 2 * np.pi
    )
    axB.ax0.loglog(F_Hz, abs(h))
    axB.ax0.loglog(F_Hz, abs(statespace.ss2xfer(
        A, B, C, D, E=E,
        F_Hz = F_Hz,
    )))
    ratio = h / statespace.ss2xfer(
        A, B, C, D, E=E,
        F_Hz = F_Hz,
    )
    axB.ax1.semilogx(F_Hz, abs(ratio))
    axB.ax1.semilogx(F_Hz, np.angle(ratio))
    axB.save(path.join(tpath, 'plot.pdf'))
    return


def test_2x2_ABCDE_c2r(plot, tpath):
    F_Hz = logspaced(.1, 100, 100)

    A = np.array([
        [1, 10],
        [0, 10]
    ])
    B = np.array([
        [1],
        [1]
    ])
    C = np.array([
        [1, 1],
    ])
    D = np.array([
        [0],
    ])
    E = np.array([
        [1j, 0],
        [0, 1j]
    ])
    A2, B2, C2, D2, E2 = SS.DSS_c2r(A, B, C, D, E, with_imag = True)
    #TODO make this an actual test
    h1 = statespace.ss2xfer(
        A, B, C, D, E=E,
        F_Hz = F_Hz,
    )
    h2 = statespace.ss2xfer(
        A2, B2, C2, D2, E=E2,
        F_Hz = F_Hz,
    )
    np.testing.assert_almost_equal(h1, h2)
    if plot:
        axB = mplfigB(Nrows = 1)
        axB.ax0.loglog(F_Hz, abs(h1))
        axB.ax0.loglog(F_Hz, abs(h2))
        axB.save(path.join(tpath, 'plot.pdf'))
    return


def test_eig_snip(tpath):
    print()
    F_Hz = logspaced(.1, 100, 100)

    Z = np.array([[0, 0], [0, 0]])
    A = np.array([[1, -1], [1, 1]])
    E = np.array([[1, 0], [0, 1]])
    A = np.block([[A, A], [-A, A]])
    E = np.block([[Z, E], [-E, Z]])

    tol = 1e-9
    v_pairs = SS.eigspaces_right(A, E, tol = tol)
    A_projections = []
    E_projections = []
    print([eig for eig, ev in v_pairs])
    for eigs, evects in v_pairs[:1]:
        #may need pivoting on to work correctly for projections
        print('eigs', eigs)
        A_project = evects[:, :1]
        #A_project = A_project / np.sum(A_project**2 , axis = 1)
        A_projections.append(A_project)

    A_projections = np.hstack(A_projections)
    Aq, Ar = scipy.linalg.qr(A_projections, mode = 'economic')
    idx_split = Aq.shape[0] - Ar.shape[0]
    A_SS_projection = np.diag(np.ones(A.shape[0])) - Aq @ Aq.T.conjugate()
    Aq, Ar, Ap = scipy.linalg.qr(A_SS_projection, mode = 'economic', pivoting = True)
    for idx in range(Aq.shape[0]):
        if np.sum(Ar[-1 - idx]**2) < tol:
            continue
        else:
            break
    idx_split = Aq.shape[0] - idx
    p_project_imU = Aq[:, :idx_split]
    p_project_kerU = Aq[:, idx_split :]

    E_projections = E @ p_project_kerU
    Eq, Er = scipy.linalg.qr(E_projections, mode = 'economic')

    E_SS_projection = np.diag(np.ones(A.shape[0])) - Eq @ Eq.T.conjugate()
    Eq, Er, Ep = scipy.linalg.qr(E_SS_projection, mode = 'economic', pivoting = True)
    p_project_im = Eq[:idx_split]
    p_project_ker = Eq[idx_split : ]

    A = p_project_im @ A @ p_project_imU
    E = p_project_im @ E @ p_project_imU

    w, vr = scipy.linalg.eig(A, E, left = False, right = True)
    print("EIGS", w)
    print("EIGV")
    print(vr)
    
    return

