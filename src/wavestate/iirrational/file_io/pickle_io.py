# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pickle


def load_pickle(fname):
    with open(fname) as F:
        fdict = pickle.load(F)
    return fdict

def write_pickle(fname, fdict):
    with open(fname, 'wb') as F:
        fdict = pickle.dump(fdict, F)
    return fdict



