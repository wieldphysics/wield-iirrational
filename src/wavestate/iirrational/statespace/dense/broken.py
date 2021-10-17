"""
"""

import numpy as np
import scipy
import scipy.signal
import scipy.linalg

from .eig_algorithms import eigspaces_right_real



def controllable_staircase(
    A, B, C, D, E,
    tol = 1e-7,
):
    """
    An attempt to do the staircase form using more general QR transforms for
    speed, not faster and very numerically aweful so far
    """
    from icecream import ic
    import tabulate
    Ninputs = B.shape[1]
    Nstates = A.shape[0]
    Nconstr = A.shape[1]
    Noutput = C.shape[0]

    print("E", tabulate.tabulate(abs(E[:9, :9])))
    BA, E = scipy.linalg.qr_multiply(
        E,
        np.hstack([B, A]),
        pivoting = False,
        mode = 'left'
    )
    print(E)
    B = BA[:, :Ninputs]
    A = BA[:, Ninputs:]
    del BA

    if False:
        #test stability
        BET, A, P = scipy.linalg.qr_multiply(
            A,
            np.hstack([B, E]).T,
            pivoting = True,
            conjugate = True,
            mode = 'right',
        )
        BE = BET.T
        B = BE[:, :Ninputs]
        E = BE[:, Ninputs:]
        E = E[:, P]
        C = C[:, P]
        del BE
        return A, B, C, D, E

    if True:
        #doesn't seem to need pivoting (which is good since it can't have it for BA)
        EAT, BA  = scipy.linalg.qr_multiply(
            np.hstack([B, A]),
            np.hstack([E, A]).T,
            pivoting = False,
            conjugate = False,
            mode = 'right'
        )
        EA = EAT.T
        B = BA[:, :Ninputs]
        E = EA[:, :Nstates]
        #A = EA[:, Nstates:]
        A = BA[:, Ninputs:]
        #del BA, EAT, EA
        #return A, B, C, D, E
        #print(tabulate.tabulate(abs(Ast[:9, :9])))
        #return A, B, C, D, E

    #Ast = Ast[::-1, ::-1]
    print()
    #A = A[::-1, ::-1]
    #E = E[::-1, ::-1]
    #B = B[::-1, :]
    #C = C[:, ::-1]
    #return A, B, C, D, E
    print(Ninputs)
    print(tabulate.tabulate(abs(E[:9, :9])))
    #A = A[::-1, ::-1]
    #E = E[::-1, ::-1]
    #B = B[::-1, :]
    #C = C[:, ::-1]
    Q, ET, P = scipy.linalg.qr(
        E.T,
        pivoting = True,
    )
    E = E[P, :]
    A = A[P, :]
    B = B[P, :]
    E = E @ Q.T
    A = A @ Q.T
    C = C @ Q.T
    #print("Q", tabulate.tabulate(abs(Q[-9:, -9:])))
    #print("A", tabulate.tabulate(abs(A[:9, :9])))
    print("E", tabulate.tabulate(abs(E[:9, :9])))
    return A, B, C, D, E

    #print("SHAPE", A.shape, C.shape)
    #ACT, ET  = scipy.linalg.qr_multiply(
    #    E.T,
    #    np.hstack([A.T, C.T]),
    #    pivoting = False,
    #    conjugate = True,
    #    mode = 'left'
    #)
    #print("1", E.shape, ET.shape)
    #E = ET.T
    #AC = ACT.T
    #print("SHAPE", AC.shape, A.shape, C.shape)
    #A = AC[:Nconstr, :]
    #C = AC[Nconstr:, :]
    #del ET, ACT, AC
    #return A, B, C, D, E
