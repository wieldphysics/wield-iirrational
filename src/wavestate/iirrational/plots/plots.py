#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import matplotlib.pyplot as plt
from wavestate import declarative

from . import plot_fit
from . import plot_residuals
from . import overlays
from . import ZP
from . import flags
from . import rel_comparison


class IIRPlots(declarative.OverridableObject):
    #this magic causes the OverridableObject constructor arguments
    #to be SAVED into _overridable_object_save_kwargs
    #this allows convenient copy/update methods
    _overridable_object_save_kwargs = True
    _overridable_object_kwargs = None

    zeroID_color = 'darkorange'
    zeroOD_color = 'darkgreen'
    poleID_color = 'blue'
    poleOD_color = 'red'
    label_data   = 'data'
    label_fit    = 'fit ({nP}P, {nZ}Z, Rsq ={Rsq:.2e})'

    show_immediate = False
    show_saved     = False
    force_save     = False

    def new(self, fname = None, **kwargs):
        kw = dict(self._overridable_object_kwargs)
        kw.update(kwargs)
        return self.__class__(**kw)

    def show(self):
        plt.show()

    def _post_plot(self, axB, fname, kwargs):
        do_show = True

        #don't show if in a subcommand (where fig is already specified to render within)
        if kwargs.get('ax_bunch', None) is not None or kwargs.get('fig', None) is not None:
            do_show = False

        if fname is not None:
            axB.save(fname)
            if not self.show_saved:
                plt.close(axB.fig)
            elif do_show and self.show_immediate:
                plt.show()
        elif self.force_save:
            raise RuntimeError("Filename to save not specified and force_save active")
        elif do_show and self.show_immediate:
            plt.show()

    def plot_fit(self, fitter, fname = None, **kwargs):
        axB = plot_fit.plot_fit(self, fitter, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_residuals(self, fitter = None, fname = None, **kwargs):
        axB = plot_residuals.plot_residuals(self, fitter, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_rel_comparison(self, fitter, fitter_ref, fname = None, **kwargs):
        axB = rel_comparison.plot_rel_comparison(self, fitter, fitter_ref, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_ZP(self, fitter, fname = None, **kwargs):
        axB = ZP.plot_ZP(self, fitter, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_ZP_S(self, fitter, fname = None, **kwargs):
        axB = ZP.plot_ZP_S(self, fitter, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_fitter_flag(self, fitter, fname = None, **kwargs):
        axB = flags.plot_fitter_flag(self, fitter, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_fitter_flag_residuals(self, fitter = None, fname = None, **kwargs):
        axB = flags.plot_fitter_flag_residuals(self, fitter = fitter, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_fitter_flag_compare(self, fitter, fitter_ref, fname = None, **kwargs):
        axB = flags.plot_fitter_flag_compare(self, fitter, fitter_ref, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB

    def plot_ZP_grab(self, fitter, duals, fname = None, **kwargs):
        axB = overlays.plot_ZP_grab(self, fitter, duals, **kwargs)
        self._post_plot(axB, fname, kwargs)
        return axB
