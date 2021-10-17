"""
"""



from .TF import (
    TF_ZPK,
    ZtoS,
    StoZ,
    SorZtoSorZ,
    ZPK_fill,
)

from .roots_bin import (
    roots_bin_type,
    roots_re_pair,
)

def abs_sq(x):
    return x.real**2 + x.imag**2

def norm_sq(v):
    return v.real.dot(v.real) + v.imag.dot(v.imag)

from numpy.polynomial.polynomial import polyvalfromroots

from .roots_matching import (
    nearest_idx,
    nearest_unique_idx,
    nearest_pairs,
    nearest_unique_pairs,
    SOS_pair_rolloff,
    match_SOS_pairs,
)

