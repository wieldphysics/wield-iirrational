#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import IPython

_ip = IPython.get_ipython()

if _ip is not None:
    _ip.magic("load_ext autoreload")
    _ip.magic("autoreload 2")

    # if this is run from the console then inline can't be found. This hack seems to get around it
    try:
        import ipykernel.pylab.backend_inline

        backend = ipykernel.pylab.backend_inline.InlineBackend.instance()
        backend.rc.clear()

        _ip.magic("matplotlib inline")
        _ip.magic("pylab inline")
    except Exception:
        _ip.magic("matplotlib")
        _ip.magic("pylab")
    # mpl.use('GTK3Cairo')

import numpy as np
import matplotlib as mpl
import matplotlib
from matplotlib import gridspec
import matplotlib.pyplot as plt

from IPython.display import (
    display,
    display_pretty,
    display_html,
    display_jpeg,
    display_png,
    display_json,
    display_latex,
    display_svg,
    display_javascript,
    display_markdown,
    FileLink,
    Image,
    SVG,
    clear_output,
    Audio,
    Javascript,
    Markdown,
)

# for more options in mpl
import wavestate.iirrational.utilities.mpl

from wavestate.iirrational.utilities.mpl.utils import (
    indexed_cmap,
)

from wavestate.iirrational.utilities.mpl import (
    AutoPlotSaver,
    mplfigB,
    asavefig,
    generate_stacked_plot_ax,
    style_6p5in,
    style_3p5in,
    style_tex_serif,
    style_serif,
    setup_log_xticks,
)


try:
    import tabulate
except ImportError:
    pass


def setup_auto_savefig(ipynb_name, check_warn=False):
    from os import path

    new = path.splitext(ipynb_name)[0] + "-ipynb"
    asavefig.org_subfolder = new


if _ip is not None:
    asavefig.fixname = True
