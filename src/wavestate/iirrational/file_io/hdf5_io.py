# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import h5py
from declarative.bunch.hdf_deep_bunch import HDFDeepBunch


def load_hdf5(fname):
    #with h5py.File(fname) as h5F:
    h5F = h5py.File(fname, 'r')
    fdict = HDFDeepBunch(h5F, writeable = False)
    return fdict

def write_hdf5(fname, fdict):
    with h5py.File(fname, 'w') as h5F:
        hdf = HDFDeepBunch(h5F, writeable = True)
        hdf.update_recursive(fdict)
    return

