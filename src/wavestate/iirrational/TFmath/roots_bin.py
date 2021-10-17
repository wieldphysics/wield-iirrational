"""
"""

import numpy as np

from . import roots_matching
conj_tol = 1e-12
real_tol = 1e-12
conj_tol = 1e-2
real_tol = 1e-4
conj_tol = 1e-6
real_tol = 1e-6

def are_conjugates(r1, r2):
    r2 = np.conjugate(r2)
    if abs(r1) < .8:
        return abs(r1 - r2) < conj_tol
    else:
        return abs((r1/r2) - 1) < conj_tol

def is_real(r1):
    if abs(r1.real) < .8:
        return abs(np.imag(r1)) < real_tol
    else:
        return abs(np.imag(r1) / np.real(r1)) < real_tol


def roots_re_pair(r_r, r_c):
    roots = list(r_r)
    for r in r_c:
        roots.append(r)
        roots.append(r.conjugate())
    return np.asarray(roots)


def roots_bin_type(
        roots,
        policy = 'auto',
        F_nyquist_Hz  = None,
        strict        = True,
        simple_output = False,
        real_tol      = 1e-6,
        conj_tol      = 1e-4,
):
    """
    roots_r are the real roots
    roots_c are just the positive complex roots
    roots_u are any unsorted roots
    """
    assert(policy in ['auto', 'pos', 'strict', 'drop_pos', 'drop_neg'])

    def are_same(r1, r2):
        if abs(r1) < .8:
            return abs(r1 - r2) < conj_tol
        else:
            return abs((r1/r2) - 1) < conj_tol

    @np.vectorize
    def are_real(r1):
        if abs(r1.real) < .8:
            return abs(np.imag(r1)) < real_tol
        else:
            return abs(np.imag(r1) / np.real(r1)) < real_tol

    roots_r = []
    roots_u = []
    seen_cplx_pos = False
    seen_cplx_neg = False

    roots_u = []
    for root in roots:
        if are_real(root):
            roots_r.append(np.real(root))
        elif np.imag(root) > 0:
            seen_cplx_pos |= True
            roots_u.append(root)
        else:
            seen_cplx_neg |= True
            roots_u.append(root)

    if policy == 'auto':
        if seen_cplx_neg and seen_cplx_pos:
            policy = 'strict'
        else:
            policy = 'pos'

    roots_c = []
    if policy == 'strict':
        roots_u = np.asarray(roots_u)
        pos_select = roots_u.imag > 0
        roots_neg = roots_u[~pos_select]
        roots_pos = roots_u[pos_select]
        rB = roots_matching.nearest_unique_pairs(roots_pos, roots_neg.conjugate())
        roots_u = list(rB.l1_remain) + [r.conjugate() for r in rB.l2_remain]
        for r1, r2 in rB.r12_list:
            if not strict or are_same(r1, r2):
                #roots_c.append((r1 + r2) / 2)
                #TODO, this seems to work better, not clear why..
                roots_c.append(r1)
            else:
                roots_u.append(r1)
                roots_u.append(r2.conjugate())
    elif policy == 'pos':
        for root in roots_u:
            if np.imag(root) > 0:
                roots_c.append(root)
            else:
                roots_c.append(root.conjugate())
            roots_u = []
    elif policy == 'drop_neg':
        for root in roots_u:
            if np.imag(root) > 0:
                roots_c.append(root)
            roots_u = []
    elif policy == 'drop_pos':
        for root in roots_u:
            if np.imag(root) < 0:
                roots_c.append(root.conjugate())
            roots_u = []

    if not simple_output:
        return (
            np.asarray(roots_r),
            np.asarray(roots_c),
            np.asarray(roots_u),
            policy
        )
    else:
        if roots_u:
            raise RuntimeError("Unbalanced or unpairable roots!")
        return roots_r, roots_c
