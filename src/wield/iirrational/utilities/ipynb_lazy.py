#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import numpy as np
import scipy
import scipy.linalg
import scipy.signal

from wield import declarative
from wield.bunch import DeepBunch
from wield.bunch.hdf_deep_bunch import HDFDeepBunch

from ..plots import (
    plot_fit,
    plot_ZP,
    plot_fitter_flag,
    plot_fitter_flag_compare,
    plot_ZP_grab,
)

from ..utilities.ipynb.displays import *

from ..utilities.np import logspaced
from ..utilities.mpl import mplfigB, generate_stacked_plot_ax

from ..fitters_ZPK import ZPKrep2MRF, MRF
from ..fitters_rational import RationalDiscFilter, ChebychevFilter

from .. import v1
from .. import v2

from ..testing import IIRrational_data


# run version printer from function to not further pollute namespace
def print_version():
    from .. import auto_version
    from .. import version

    print(
        "wield.iirrational version: {} (git:{})".format(
            version, auto_version.git_shorthash
        )
    )


print_version()
del print_version
