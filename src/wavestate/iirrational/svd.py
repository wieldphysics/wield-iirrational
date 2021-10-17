# -*- coding: utf-8 -*-
"""
"""


import scipy
import scipy.linalg

try:
    #if the package was built then we can go faster
    from .lapack.svd import svd as svd_fast
except ImportError:
    svd_fast = None


def SVD_SV(
    a,
    overwrite_a = False,
    n_smallest  = None,
):
    if svd_fast is not None:
        if n_smallest is not None:
            S, V = svd_fast(
                a,
                compute_u   = False,
                compute_v   = True,
                il          = -n_smallest,
                overwrite_a = overwrite_a,
            )
        else:
            S, V = svd_fast(
                a,
                compute_u   = False,
                compute_v   = True,
                overwrite_a = overwrite_a,
            )
    else:
        U, S, V = scipy.linalg.svd(
            a,
            compute_uv = True,
            full_matrices = False,
            overwrite_a = overwrite_a,
        )
    return S, V

