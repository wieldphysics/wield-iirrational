#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

#import warnings
import contextlib
from wavestate import declarative

from .. import fitters_ZPK
from .. import annotate
from .. import plots
from . import algorithms


def fitter_order(fitter):
    return max(len(fitter.zeros), len(fitter.poles))


class FitterUnacceptable(Exception):
    pass


class FitAid(object):

    def __init__(
        self,
        doc_db,
        hints = None,
        residuals_type = 'log',
    ):
        self.fitter = None
        self.fitter_current    = None
        self.fitter_lowres_avg = None
        self.fitter_lowres_max = None
        self.fitter_lowres_med = None
        self.fitter_loword     = None

        self.doc_db_method     = None

        if hints is None:
            hints = dict()

        if isinstance(hints, (list, tuple)):
            usehints = dict()
            for hdict in hints:
                usehints.update(hdict)
        else:
            usehints = dict(hints)

        self.hints = usehints
        self.hints_defaults = dict()

        self.doc_db = doc_db
        self.doc_stack = []

        self.section = [0]
        self.section_stack = []

        #TODO make this a hint
        self.verbosity = 5

        return

    @property
    def residuals_type(self):
        return self.hint('residuals_type', default = 'log')

    def _check_residuals(self, fitter):
        if not isinstance(fitter, fitters_ZPK.MultiReprFilterBase):
            raise FitterUnacceptable()
        if fitter.residuals_type == self.residuals_type:
            if fitter.residuals_log_im_scale != 1:
                f_new = fitter.copy(
                    residuals_type = self.residuals_type,
                    residuals_log_im_scale = 1,
                )
                return f_new
        return None

    def hint_set(self, hname, hval, default = False):
        if default:
            self.hints.setdefault(hname, hval)
        else:
            self.hints[hname] = hval
        return

    def hint_arg(self, func_arg, *args, **kwargs):
        if func_arg is None:
            return self.hint(*args, **kwargs)
        else:
            return func_arg

    def hint(self, *args, **kwargs):
        superarg = []
        for arg in args:
            if isinstance(arg, (list, tuple)):
                superarg.extend(arg)
            else:
                superarg.append(arg)

        default = kwargs.get('default', declarative.NOARG)
        if default is declarative.NOARG:
            raise RuntimeError("Must have default return")

        first_hint = superarg[0]

        lh_set = self.hints_defaults.get(first_hint, None)
        if lh_set is None:
            lh_set = set()
            self.hints_defaults[first_hint] = lh_set
        lh_set.add(default)

        for key in superarg:
            key = key.format(**kwargs)
            ret = self.hints.setdefault(key, None)
            if ret is not None:
                return ret
        return default

    def fitter_check(
        self,
        fitter_new = None,
        hint_name = None,
        update    = True,
        variant   = None,
        validate  = True,
        annotate  = "unknown",
    ):
        """
        """

        if variant is None:
            raise RuntimeError('Must specify variant for fitter_check, such as ["OrdUp", "OrdDn", "OrdC"]')

        if self.fitter is None:
            if update:
                self.fitter = fitter_new
                self.fitter_update()
            else:
                return True

        if fitter_new is None:
            fitter_new = self.fitter

        #TODO, should check that it is using standard residuals

        resavgN = lambda : fitter_new.residuals_average
        resmaxN = lambda : fitter_new.residuals_max
        resmedN = lambda : fitter_new.residuals_med
        resavgC = lambda : self.fitter_current.residuals_average
        resmaxC = lambda : self.fitter_current.residuals_max
        resmedC = lambda : self.fitter_current.residuals_med
        resavgB = lambda : self.fitter_lowres_avg.residuals_average
        resmaxB = lambda : self.fitter_lowres_max.residuals_max
        resmedB = lambda : self.fitter_lowres_med.residuals_med

        print('resavgC:', resavgC())

        if hint_name is None:
            hint_name = 'default'

        improved = [True]

        def check_hint(
            hints,
            calc,
        ):
            hintval = self.hint(hints, default = None)
            if hintval is not None and calc > hintval:
                self.log("Failed {}, {} > {}".format(hints[0], calc, hintval))
                improved[0] = False

        check_hint([
            'resavg_Ethresh{0}.{1}'.format(variant, hint_name),
            'resavg_Ethresh{0}'.format(variant),
        ], resavgN())

        check_hint([
            'resmax_Ethresh{0}.{1}'.format(variant, hint_name),
            'resmax_Ethresh{0}'.format(variant),
        ], resmaxN())

        check_hint([
            'resavg_Rthresh{0}.{1}'.format(variant, hint_name),
            'resavg_Rthresh{0}'.format(variant),
        ], (resavgN() / resavgC()))

        check_hint([
            'resmax_Rthresh{0}.{1}'.format(variant, hint_name),
            'resmax_Rthresh{0}'.format(variant),
        ], (resmaxN() / resmaxC()))

        check_hint([
            'resmed_Rthresh{0}.{1}'.format(variant, hint_name),
            'resmed_Rthresh{0}'.format(variant),
        ], (resmedN() / resmedC()))

        check_hint([
            'resavgB_Rthresh{0}.{1}'.format(variant, hint_name),
            'resavgB_Rthresh{0}'.format(variant),
        ], (resavgN() / resavgB()))

        check_hint([
            'resmaxB_Rthresh{0}.{1}'.format(variant, hint_name),
            'resmaxB_Rthresh{0}'.format(variant),
        ], (resmaxN() / resmaxB()))

        check_hint([
            'resmedB_Rthresh{0}.{1}'.format(variant, hint_name),
            'resmedB_Rthresh{0}'.format(variant),
        ], (resmedN() / resmedB()))

        self.log('IMPROVED: ', improved[0])
        if improved[0] and update:
            self.fitter_update(
                fitter_new,
                annotate = None,
            )
            if annotate is not None:
                self.annotate(
                    name      = 'fitter_check (improved)',
                    about     = annotate,
                    fitter    = fitter_new,
                    plotter   = plots.plot_fitter_flag,
                    method    = self.doc_db_method,
                    verbosity = 6,
                )
        else:
            if annotate is not None:
                self.annotate(
                    name      = 'fitter_check (not improved)',
                    fitter    = fitter_new,
                    about     = annotate,
                    plotter   = plots.plot_fitter_flag,
                    method    = self.doc_db_method,
                    verbosity = 10,
                )
        val_func = self.hint('fitter_check_validate', default = None)
        if validate and val_func is not None and fitter_new is not None:
            val_func(self, fitter_new)
        return improved[0]

    def fitter_checkup(self):
        improved = self.fitter_check(
            self.fitter,
            variant = 'OrdC',
        )
        if not improved:
            self.fitter = self.fitter_current.copy()
        return improved

    def fitter_update(
            self,
            fitter = None,
            validate = True,
            annotate  = "unknown",
    ):
        if fitter is None:
            fitter = self.fitter
        else:
            self.fitter = fitter
        avg1 = fitter.residuals_average
        algorithms.sign_check_flip(fitter)
        avg2 = fitter.residuals_average
        assert(avg2 <= avg1)

        if annotate is not None:
            self.annotate(
                name      = 'fitter_update',
                about     = annotate,
                fitter    = fitter,
                plotter   = plots.plot_fitter_flag,
                method    = self.doc_db_method,
                verbosity = 8,
            )

        val_func = self.hint('fitter_update_validate', default = None)
        if validate and val_func is not None and fitter is not None:
            val_func(self, fitter)

        try:
            fitter_use = self._check_residuals(fitter)
        except FitterUnacceptable:
            return

        if fitter_use is not None:
            fitter_use.optimize(aid = self)
            fitter_copy = [fitter_use]
        else:
            fitter_use = fitter
            fitter_copy = []

        def copy1():
            if not fitter_copy:
                fitter_copy.append(fitter.copy())
            return fitter_copy[0]

        #print("ZPK", fitter_use.ZPKsf)

        new_res = fitter_use.residuals_average
        new_res_max = fitter_use.residuals_max
        new_res_med = fitter_use.residuals_med
        new_ord = fitter_order(fitter_use)

        if self.fitter_current is not None:
            old_res     = self.fitter_lowres_avg.residuals_average
            old_res_max = self.fitter_lowres_max.residuals_max
            old_res_med = self.fitter_lowres_med.residuals_med
            old_ord     = fitter_order(self.fitter_loword)

            if new_res < old_res:
                self.fitter_lowres_avg = copy1()

            if new_res_max < old_res_max:
                self.fitter_lowres_max = copy1()

            if new_res_med < old_res_med:
                self.fitter_lowres_med = copy1()

            if new_ord < old_ord:
                self.fitter_loword = copy1()

            self.fitter_current = copy1()
        else:
            self.fitter_current    = copy1()
            self.fitter_lowres_avg = copy1()
            self.fitter_lowres_max = copy1()
            self.fitter_lowres_med = copy1()
            self.fitter_loword     = copy1()

        #print("UPDATED: ", self.fitter_current.residuals_average)
        
        return fitter_use

    def log(self, *args, **kwargs):
        #hint = kwargs.get('hint', None)

        if self.hint('no_log', default = False):
            return

        if self.hint('log_print', default = True):
            print(*args)

        #TODO put in the doc_db

    def digest_write(
        self,
        folder,
        format = 'markdown',
        **kwargs
    ):
        if format in ['md', 'markdown']:
            return annotate.docdb2markdown(
                folder,
                self.doc_db,
                **kwargs
            )
        else:
            raise RuntimeError("Format Not Recognized, must be one of: ['md', 'markdown']")

    @contextlib.contextmanager
    def annotate_into(self, name):
        doc_db_prev = self.doc_db

        self.doc_db = annotate.annotate_into(
            self.doc_db,
            name = name,
        )
        yield

        self.doc_db = doc_db_prev
        return

    @contextlib.contextmanager
    def annotate_method(self, name):
        doc_db_method_prev = self.doc_db_method

        self.doc_db_method = name
        yield

        self.doc_db_method = doc_db_method_prev
        return

    def annotate(
        self,
        name,
        verbosity,
        about,
        method          = None,
        verbosity_limit = None,
        fitter          = None,
        plotter         = None,
    ):
        if method is None:
            method = self.doc_db_method
        annotate.annotate(
            self.doc_db,
            name            = name,
            verbosity       = verbosity,
            method          = method,
            about           = about,
            verbosity_limit = verbosity_limit,
            fitter          = fitter,
            plotter         = plotter,
        )

