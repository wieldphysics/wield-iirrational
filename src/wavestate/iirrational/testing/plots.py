"""
Some utilities to help with dataset creation and annotation
"""

import os
import os.path as path
import contextlib

from IIRrational import (
    fitters_rational,
    fitters_ZPK,
    representations,
)

#import matplotlib.pyplot as plt
from IIRrational import plots
from IIRrational.annotate.exporters import npm_markdown_pdf

@contextlib.contextmanager
def plot_on_assert(testfile, request, fitter, plot_anyway = False):
    if isinstance(request, (str)):
        rname = request
    else:
        rname = request.node.name
    fname = path.join(
        path.split(testfile)[0],
        'plots/{}.pdf'.format(rname)
    )

    def plot():
        if isinstance(
                fitter,
                (
                    fitters_rational.DataFilterBase,
                    fitters_ZPK.MultiReprFilterBase,
                    representations.ZPKwData,
                )
        ):
            axB = plots.plot_fitter_flag(
                fitter,
                fname = fname
            )
        else:
            axB = fitter()
            axB.save(fname)
    try:
        yield
    except AssertionError:
        plot()
        raise

    if plot_anyway:
        plot()
    else:
        try:
            os.remove(fname)
        except OSError:
            pass
    return


@contextlib.contextmanager
def digest_on_assert(testfile, request, aid, plot_anyway = False):
    if isinstance(request, (str)):
        rname = request
    else:
        rname = request.node.name
    fname = path.join(
        path.split(testfile)[0],
        'plots/{}/'.format(rname)
    )

    def plot():
        aid.digest_write(
            fname,
            plot_verbosity     = 9,
            clear_plots        = True,
            plot_format        = 'png',
            dpi                = 200,
            exporter           = npm_markdown_pdf('digest.pdf'),
            MP_workers         = 1,
        )
    try:
        yield
    except AssertionError:
        plot()
        raise

    if plot_anyway:
        plot()
    else:
        try:
            #os.remove(fname)
            pass
        except OSError:
            pass
    return
