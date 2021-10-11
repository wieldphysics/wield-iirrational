# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
from os import path


def relpy(pyfile, fpath):
    return path.join(path.split(pyfile)[0], fpath)


from .fake_aid import ensure_aid
