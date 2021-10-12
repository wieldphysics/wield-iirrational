# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import declarative
import os
import re


ext2type = {
    '.h5'      : 'hdf5',
    '.hdf'     : 'hdf5',
    '.hdf5'    : 'hdf5',
    '.pkl'     : 'pickle',
    '.pickle'  : 'pickle',
    '.yaml'    : 'yaml',
    '.yml'     : 'yaml',
    '.json'    : 'json',
    '.ini'     : 'ini',
    '.mat'     : 'mat',
    '.m'       : 'mat',
    '.txt'     : 'csv',
    '.csv'     : 'csv',
    '.txt.gz'  : 'csv',
    '.csv.bz2' : 'csv',
}


type2type = {
    'h5'     : 'hdf5',
    'hdf'    : 'hdf5',
    'hdf5'   : 'hdf5',
    'yaml'   : 'yaml',
    'yml'    : 'yaml',
    'json'   : 'json',
    'ini'    : 'ini',
    'matlab' : 'mat',
    'mat'    : 'mat',
    'm'      : 'mat',
    'csv'    : 'csv',
    'pkl'    : 'pickle',
    'pickle' : 'pickle',
    'special' : 'special',
}


type2features = {
    'hdf5'   : dict(
        data    = True,
        config  = True,
        compact = True,
        complex = True,
        ndarray = True,
    ),
    'yaml'   : dict(
        data    = True,
        config  = True,
        compact = False,
        complex = False,
        ndarray = False,
    ),
    'json'   : dict(
        data    = True,
        config  = False,
        compact = False,
        complex = False,
        ndarray = False,
    ),
    'ini'    : dict(
        data    = False,
        config  = True,
        compact = False,
        complex = False,
        ndarray = False,
    ),
    'mat'    : dict(
        data    = True,
        config  = True,
        compact = True,
        complex = True,
        ndarray = True,
    ),
    'csv'    : dict(
        data    = True,
        config  = False,
        compact = True,
        complex = False,
        ndarray = True,
    ),
    'pickle' : dict(
        data    = True,
        config  = True,
        compact = True,
        complex = True,
        ndarray = True,
    ),
}

re_FILEKEY = re.compile(r'(.*)\[(.*)\]')


def determine_type(fname):
    fspl = fname.split('::')

    if len(fspl) == 1:
        fname, = fspl
        m = re_FILEKEY.match(fname)
        if m:
            fname = m.group(1)
            subkey = m.group(2)
        else:
            subkey = None
        fbase, fext = os.path.splitext(fname)
        if fname[0] == ':' and fname[-1] == ':':
            ftype = 'special'
            fname = fname[1:-1]
        else:
            ftype = ext2type[fext.lower()]
    elif len(fspl) == 2:
        fname, ftype = fspl
        if m:
            fname = m.group(1)
            subkey = m.group(2)
        else:
            subkey = None
        ftype = type2type[ftype.lower()]
        if (
                ftype == 'special'
                and fname[0] == ':'
                and fname[-1] == ':'
        ):
            fname = fname[1:-1]
    else:
        #TODO, be more specific about which argument caused this
        raise RuntimeError("Too many '::' in filetype specification")

    if ftype != 'special':
        fname = os.path.abspath(fname)

    return wavestate.bunch.Bunch(
        fname  = fname,
        subkey = subkey,
        ftype  = ftype,
    )


