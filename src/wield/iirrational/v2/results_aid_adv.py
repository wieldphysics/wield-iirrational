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
from wield.bunch import Bunch

# TODO, package with lib


# from .arguments import ArgumentError
# from .. import fitters_ZPK

from ..external.tabulate import tabulate
from .. import plots
from . import results_aid


class ResultsAidAdv(results_aid.ResultsAid):
    def __init__(self, fit_aid, **kwargs):
        super(ResultsAidAdv, self).__init__(fit_aid, **kwargs)
        self._investigations = {}
        self._investigations.update(fit_aid.investigations)
        for name, func in self._investigations.items():
            setattr(
                self,
                name,
                (lambda func: lambda *args, **kwargs: func(self, *args, **kwargs))(
                    func
                ),
            )
        return

    def investigate_fit_plot(self, fname=None, xscale="log_zoom", **kwargs):
        """
        Generates a plot of the fit. Returns a Bunch object which
        stores the axes and figure of the plot.
        """
        # axB = plots.plot_fitter_flag(fitter=self.fitter, xscale=xscale, **kwargs)
        axB = plots.plot_fitter_flag_residuals(fitter=self.fitter, xscale=xscale, **kwargs)
        if fname is not None:
            axB.save(fname)
        return axB

    def investigate_order_arrays(self, local_best=True):
        obo, rbo, fbo, lbo = self._fitters_by_order

        # only returns or displays those which are the best for that order or below
        if local_best:
            obo = obo[lbo]
            rbo = rbo[lbo]
            fbo = [fbo[idx] for idx in lbo]

        rmaxbo = []
        rmedbo = []
        for fit in fbo:
            rmaxbo.append(fit.residuals_max)
            rmedbo.append(fit.residuals_med)
        return Bunch(
            order=obo,
            res_avg=rbo,
            res_max=rmaxbo,
            res_med=rmedbo,
        )

    def investigate_order_console(
        self,
        tablefmt="simple",
        print_function=print,
    ):
        oB = self.investigate_order_arrays()
        pval = tabulate(
            np.asarray([oB.order, oB.res_avg, oB.res_med, oB.res_max]).T,
            ["max(z, p)\norder", "ChiSq.\navg. res.", "\nmed. res.", "\nmax. res."],
            tablefmt=tablefmt,
        )
        if print_function is not None:
            print_function(pval)
        else:
            return pval

    def investigate_order_plot(self, fname=None, xscale="log_zoom", **kwargs):
        obo, rbo, fbo, lbo = self._fitters_by_order
        axB = plots.plot_fitter_flag_residuals(
            fitter=self.fitter, xscale=xscale, **kwargs
        )
        for idx in reversed(range(len(fbo))):
            fit = fbo[idx]
            if fit is self.fitter:
                continue
            if not lbo[idx]:
                continue
            axB.plot_fit(
                fit,
                label=("order: {}, ChiSq: {:.3e}").format(
                    fit.order, fit.residuals_average
                ),
            )
        if fname is not None:
            axB.save(fname)
        else:
            axB.finalize()
        return axB

    def investigate_covarianceRI(self, output="print"):
        """
        outputs the covariance matrix for the real/imaginary parts of the chosen
        filter.
        """
        # TODO
        raise NotImplementedError()

    def investigate_valuesRI_print(self, output="print"):
        """
        outputs the covariance matrix for the real/imaginary parts of the chosen
        filter.

        output specifies the output format, which may be one of:
        - "print":
        - "arguments":
        - "matrix":
        - "values":
        """
        # TODO
        raise NotImplementedError()

    def investigate_fisherRI_print(self, output="print"):
        """
        outputs the covariance matrix for the real/imaginary parts of the chosen
        filter.

        output specifies the output format, which may be one of:
        - "print":
        - "arguments":
        - "matrix":
        - "values":
        """
        # TODO
        raise NotImplementedError()

    # def digest(
    #    self,
    #    fname = None,
    #    folder = None,
    #    verbosity = None,
    #    **kwargs
    # ):
    #    """
    #    Generates a document with the logging output, intermediate fits and
    #    investigations related to the fitting call. This is verbose and slow,
    #    but useful.

    #    The default is to print and plot to the output,
    #    if running from an ipython/jupyter notebook.
    #    """
    #    raise NotImplementedError()
