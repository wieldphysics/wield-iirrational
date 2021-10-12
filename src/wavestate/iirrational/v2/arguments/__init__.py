# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import declarative

from .base import (
    ArgumentError,
    grab_kwargs,
    grab_kwarg_hints,
    kw_ZPKrep_build,
    check_remaining_arguments,
    transfer_kw,
)
from . import logging
from . import standardargs
from . import adjustments
from . import ranges

kw_hints = wavestate.bunch.Bunch()
kw_hints.update(logging.kw_hints)
kw_hints.update(standardargs.kw_hints)
kw_hints.update(adjustments.kw_hints)
kw_hints.update(ranges.kw_hints)

__all__ = [
    ArgumentError,
    grab_kwargs,
    grab_kwarg_hints,
    kw_ZPKrep_build,
    kw_hints,
    logging,
    standardargs,
    adjustments,
    ranges,
]
