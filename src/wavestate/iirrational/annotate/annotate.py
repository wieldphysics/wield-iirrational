# -*- coding: utf-8 -*-
"""
"""

from wavestate.bunch import DeepBunch
import os
import re
import errno
import os.path as path
import time

def create(
        doc_db,
        name,
        verbosity_limit = 3,
):
    if doc_db is None:
        doc_db = DeepBunch()
        doc_db.verbosity_limit = verbosity_limit
    doc_db.name = name
    doc_db.time = time.time()
    doc_db.section = []
    return doc_db


def annotate(
    doc_db,
    name,
    verbosity,
    method,
    about,
    verbosity_limit = None,
    fitter          = None,
    plotter         = None,
):
    if name is None:
        doc = doc_db
        if doc is not None:
            doc.method    = method
            doc.about     = about
            doc.verbosity = verbosity
            doc.time = time.time()
            if not doc.verbosity:
                doc.verbosity = verbosity
            else:
                verbosity = doc.verbosity
            if doc.verbosity_limit:
                verbosity_limit = doc.verbosity_limit
            else:
                verbosity_limit = 10
            doc.plotter = plotter
            if verbosity_limit >= verbosity:
                if fitter is not None:
                    doc.fitter = fitter.copy()
                return doc
    else:
        if doc_db is not None:
            doc = doc_db[name]
            doc.name      = name
            doc.method    = method
            doc.about     = about
            doc.time = time.time()
            if not doc.verbosity:
                doc.verbosity = verbosity
            else:
                verbosity = doc.verbosity
            if verbosity_limit is None:
                verbosity_limit = doc.verbosity_limit
            if not verbosity_limit:
                verbosity_limit = 10
            if not doc.verbosity_limit:
                doc.verbosity_limit = verbosity_limit
            else:
                verbosity_limit = doc.verbosity_limit
            doc.plotter = plotter
            if verbosity_limit >= verbosity:
                seq = doc_db.sequence
                if not seq:
                    seq = []
                    doc_db.sequence = seq
                doc.section = doc_db.section + [len(seq)]
                seq.append(name)
                if fitter is not None:
                    doc.fitter = fitter.__class__(copy = fitter)
                return doc


def annotate_into(
    doc_db,
    name,
):
    if doc_db is not None:
        seq = doc_db.sequence
        if not seq:
            seq = []
            doc_db.sequence = seq
        seq.append(name)
        doc = doc_db[name]
        doc.name = name
        doc.section = doc_db.section + [len(doc_db.sequence)]
        if doc_db.verbosity_limit:
            doc.verbosity_limit = doc_db.verbosity_limit
        return doc
