#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import numpy as np
import scipy
import scipy.signal
from IIRrational import TFmath


def eigspaces_right(A, B=None, tol=1e-9):
    ###############################
    # eigenvalue subspace collection
    w, vr = scipy.linalg.eig(A, B, left=False, right=True)

    w_collect_val = [[_] for _ in w]
    w_collect_idx = [[_] for _, v in enumerate(w)]

    collected_any = True
    w = list(w)
    while collected_any:
        # print(w / (2*np.pi))
        # print(vr)
        w_near = TFmath.nearest_idx(w)
        tol = 1e-7
        collected_any = False
        for w_idx, (w_val, w_near_idx) in enumerate(zip(w, w_near)):
            if w_near_idx is None:
                continue
            w_val2 = w[w_near_idx]
            if w_val is None or w_val2 is None:
                continue
            if abs(w_val - w_val2) < tol:
                w[w_near_idx] = None
                w_collect_val[w_idx] += w_collect_val[w_near_idx]
                w_collect_val[w_near_idx] = []
                w_collect_idx[w_idx] += w_collect_idx[w_near_idx]
                w_collect_idx[w_near_idx] = []
                collected_any = True
    # print(w_collect_val)
    # print(w_collect_idx)
    w_pairs = [p for p in zip(w_collect_idx, w_collect_val) if len(p[0]) > 0]
    v_pairs = []
    # u, s, v = scipy.linalg.svd(vr)
    # print(s)
    for idxs, eigs in w_pairs:
        evects = vr[:, idxs]
        # u, s, v = scipy.linalg.svd(evects)
        # print(s)
        v_pairs.append((eigs, evects))
    # the evects output is rows are A-space, columns are eig-idx-space
    return v_pairs


def eigspaces_right_real(A, B=None, tol=1e-9):
    v_pairs = eigspaces_right(A, B=B, tol=tol)
    v_pairs_re = []
    v_pairs_im = []
    w_im = []
    for eigs, evects in v_pairs:
        eigv = np.mean(eigs)
        if abs(eigv.imag) < tol:
            # remove the imaginary part to the eigenvectors
            assert np.all(np.sum(evects.imag ** 2, axis=1) < tol)
            v_pairs_re.append((eigs, evects.real))
            continue
        v_pairs_im.append((eigs, evects))
        w_im.append(eigv)

    # this finds conjugate pairs
    w_near = TFmath.nearest_idx(w_im)
    v_pairs_im2 = []
    for idx_fr, idx_to in enumerate(w_near):
        if idx_to is None or w_near[idx_to] is None:
            continue
        if idx_fr is None or w_near[idx_fr] is None:
            continue
        if w_near[idx_to] != idx_fr:
            continue
        # unique conjugate pair
        w_near[idx_to] = None
        w_near[idx_fr] = None
        eigs1, eigv1 = v_pairs_im[idx_to]
        eigs2, eigv2 = v_pairs_im[idx_fr]
        if w_im[idx_to].imag > 0:
            eigs_use = eigs1
        else:
            eigs_use = eigs2
        v_pairs_im2.append((eigs_use, np.hstack([eigv1, eigv2])))

    v_pairs_im3 = []
    for eigs, evects in v_pairs_im2:
        # TODO, it may be the the SVD should be used here
        # this may rely on r being rank-revealing
        q, r = scipy.linalg.qr(evects.imag.T)
        # check that the imaginary space is actually reduced
        idx_cut = r.shape[0] // 2
        assert np.all(np.sum(r[idx_cut:] ** 2, axis=1) < tol)
        # same check as above (redundant)
        assert np.all(np.sum((evects.imag @ q[idx_cut:].T) ** 2, axis=1) < tol)
        # now formulate the real projection
        evects2 = evects.real @ q[idx_cut:].T

        # normalize the eigenvectors again
        evects2 = evects2 / (np.sum(evects2 ** 2, axis=0)) ** 0.5

        v_pairs_im3.append((eigs, evects2))

    return v_pairs_re + v_pairs_im3
