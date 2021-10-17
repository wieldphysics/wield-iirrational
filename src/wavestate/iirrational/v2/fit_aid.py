#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import sys
import time
import logging
import numpy as np

import contextlib
import declarative
from declarative.utilities.future_from_2 import str, unicode

from .. import fitters_ZPK
from .. import representations
from ..utilities.strings import padding_remove
from . import algorithms


class FitterUnacceptable(Exception):
    pass


class FitAid(object):
    def __init__(
        self,
        hints = None,
        hints_seen = None,
        active = True,
    ):
        self.fitter             = None
        self.fitter_special     = wavestate.bunch.Bunch()
        self._fitter_current    = None
        self._fitter_lowres_avg = None
        self._fitter_lowres_max = None
        self._fitter_lowres_med = None

        self.mtime_start = time.time()

        if hints is None:
            hints = dict()

        if isinstance(hints, (list, tuple)):
            usehints = dict()
            for hdict in hints:
                usehints.update(hdict)
        else:
            usehints = dict(hints)

        self.hints = usehints
        self.hints_seen = hints_seen

        self.section = [0]
        self.section_stack = []

        #TODO make this a hint
        self.verbosity = 5
        #self.verbosity_info
        #self.verbosity_alert
        #self.verbosity_warn
        #self.verbosity_debug
        self.active = active

        #holds a heading for the logging, as well as sets tabbing
        self.log_header_stack = ()
        #indicates how far into the header has been printed yet.
        #for the live log
        self.log_header_printed = 0

        #log entries
        self._logs = []
        self.log_number = 0

        #stores all representative fitters seen and log metadata and checkpoint metadata
        self._fitters = []
        self._fitter_factors = []

        #stores all _checkpoints and metadata
        self._checkpoints = []
        self.checkpoint_current = None

        #investigations to view
        self.investigations = dict()
        return

    def __bool__(self):
        return self.active

    def __nonzero__(self):
        return self.active

    @property
    def residuals_type(self):
        return self.hint('residuals_type', 'log')

    @property
    def residuals_type_alt(self):
        return self.hint('residuals_type_alt', 'dualA')

    def _check_residuals(self, fitter):
        """
        Get a fitter with representative residuals for comparison and storage.
        This applies all current factorizations.
        """
        if not isinstance(fitter, fitters_ZPK.MultiReprFilterBase):
            raise FitterUnacceptable()

        if self._fitter_factors:
            first_factor = self._fitter_factors[0]
            zeros = first_factor.zeros
            poles = first_factor.poles
            gain  = first_factor.gain
            for factor in self._fitter_factors[1:]:
                zeros = zeros * factor.zeros
                poles = zeros * factor.poles
                gain  = factor.gain
            xrep = representations.ZPKwData(
                ZPKrep = first_factor.ZPKrep,
                zeros  = fitter.zeros * zeros,
                poles  = fitter.poles * poles,
                gain   = fitter.gain  * gain,
            )
            #TODO, check that the representation is preserved
            return fitter.regenerate(
                ZPKrep = xrep,
                residuals_log_im_scale = 1,
            )

        #TODO, should this assert if residuals_type is wrong?
        if fitter.residuals_type == self.residuals_type:
            if np.all(fitter.residuals_log_im_scale) != 1:
                f_new = fitter.copy(
                    residuals_type = self.residuals_type,
                    residuals_log_im_scale = 1,
                )
                return f_new
        return None

    #declarative
    def fitter_orders(self, fitter = None):
        if fitter is None:
            fitter = self.fitter
        if self._fitter_factors:
            first_factor = self._fitter_factors[0]
            total_z = len(first_factor.zeros) + len(first_factor.zeros_overlay)
            total_p = len(first_factor.poles) + len(first_factor.poles_overlay)
            for factor in self._fitter_factors[1:]:
                total_z += len(factor.zeros)
                total_p += len(factor.poles)
            if fitter != 'factors':
                total_z += len(fitter.zeros)
                total_p += len(fitter.poles)
        else:
            if fitter != 'factors':
                total_z = len(fitter.zeros) + len(fitter.zeros_overlay)
                total_p = len(fitter.poles) + len(fitter.poles_overlay)
            else:
                total_z = 0
                total_p = 0
        return wavestate.bunch.Bunch(
            z      = total_z,
            p      = total_p,
            maxzp  = max(total_z, total_p),
            total  = total_z + total_p,
            reldeg = total_z - total_p,
        )

    #def hint_set(self, hname, hval, default = False):
    #    if default:
    #        self.hints.setdefault(hname, hval)
    #    else:
    #        self.hints[hname] = hval
    #    return

    def hint_has(self, hname):
        return hname in self.hints

    def hint_setdefault(self, hname, hval):
        self.hints.setdefault(hname, hval)
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

        #helper for when the known keys are being indexed
        if self.hints_seen is not None:
            keys_remapped = [key.format(**kwargs) for key in superarg]
            for idx, key in enumerate(keys_remapped):
                overrides, relateds = self.hints_seen.setdefault(key, (set(), set()))
                relateds.update(keys_remapped[:idx])
                overrides.update(keys_remapped[idx+1:])

        for key in superarg:
            key = key.format(**kwargs)
            ret = self.hints.get(key, declarative.NOARG)
            if ret is not declarative.NOARG:
                return ret
        return kwargs['default']

    def fitter_check(
        self,
        fitter_new = None,
        hint_name = None,
        update    = True,
        variant   = None,
        validate  = True,
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
        resavgC = lambda : self._fitter_current.residuals_average
        resmaxC = lambda : self._fitter_current.residuals_max
        resmedC = lambda : self._fitter_current.residuals_med
        resavgB = lambda : self._fitter_lowres_avg.residuals_average
        resmaxB = lambda : self._fitter_lowres_max.residuals_max
        resmedB = lambda : self._fitter_lowres_med.residuals_med

        if hint_name is None:
            hint_name = 'default'

        improved = [True]

        def check_hint(
            hints,
            calc,
        ):
            hintval = self.hint(hints, default = None)
            if hintval is not None and calc > hintval:
                self.log_debug(10, "Failed {}, {} > {}".format(hints[0], calc, hintval))
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

        self.log_debug(
            10,
            'IMPROVED: ', improved[0],
            'residuals: ', resavgN(), resavgC()
        )

        val_func = self.hint('fitter_check_validate', default = None)
        if validate and val_func is not None and fitter_new is not None:
            val_func(self, fitter_new)

        if improved[0] and update:
            self.fitter_update(
                fitter_new,
                representative = True,
            )
        return improved[0]

    def fitter_checkup(self):
        ord_chg = self.fitter.order_total - self._fitter_current.order_total
        if ord_chg < 0:
            variant = 'OrdDn'
        elif ord_chg > 0:
            variant = 'OrdUp'
        else:
            variant = 'OrdC'
        algorithms.sign_check_flip(self.fitter)
        improved = self.fitter_check(
            self.fitter,
            variant = variant,
        )
        if not improved:
            self.fitter = self._fitter_current.copy()
        return improved

    def fitter_update(
        self,
        fitter         = None,
        representative = False,
        validate       = True,
    ):
        if fitter is None:
            fitter = self.fitter
        else:
            self.fitter = fitter

        avg1 = fitter.residuals_average
        algorithms.sign_check_flip(fitter)
        avg2 = fitter.residuals_average
        assert(avg2 <= avg1)

        val_func = self.hint('fitter_update_validate', default = None)
        if validate and val_func is not None and fitter is not None:
            val_func(self, fitter)

        #this function will apply the de-factorizations, to create an
        #unfactored fit
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
        if representative:
            fitter_meta = wavestate.bunch.Bunch()
            fitter_meta.fitter         = copy1()
            fitter_meta.log_idx        = self.log_number
            fitter_meta.checkpoint_idx = len(self._checkpoints) - 1
            fitter_meta.valid          = True
            self._fitters.append(fitter_meta)

            new_res = fitter_use.residuals_average
            new_res_max = fitter_use.residuals_max
            new_res_med = fitter_use.residuals_med

            if self._fitter_current is not None:
                old_res     = self._fitter_lowres_avg.residuals_average
                old_res_max = self._fitter_lowres_max.residuals_max
                old_res_med = self._fitter_lowres_med.residuals_med

                if new_res < old_res:
                    self._fitter_lowres_avg = copy1()

                if new_res_max < old_res_max:
                    self._fitter_lowres_max = copy1()

                if new_res_med < old_res_med:
                    self._fitter_lowres_med = copy1()

                self._fitter_current = copy1()
            else:
                self._fitter_current    = copy1()
                self._fitter_lowres_avg = copy1()
                self._fitter_lowres_max = copy1()
                self._fitter_lowres_med = copy1()

        self.fitter_special = wavestate.bunch.Bunch()
        return fitter_use

    def factorization_push(
        self,
        fitter_factor = None,
        data_mod      = False,
    ):
        if fitter_factor is None:
            fitter_factor = self.fitter
        self._fitter_factors.append(fitter_factor)

        if data_mod:
            xrep = representations.ZPKwData(
                ZPKrep = fitter_factor.ZPKrep,
                data = self.fitter.data / fitter_factor.xfer_fit,
                gain  = 1,
                poles = (),
                zeros = (),
                poles_overlay = (),
                zeros_overlay = (),
                delay_s = self.fitter.delay_s - fitter_factor.delay_s,
            )
            self.fitter = self.fitter.regenerate(ZPKrep = xrep)
        else:
            xrep = representations.ZPKwData(
                ZPKrep = fitter_factor.ZPKrep,
                gain  = 1,
                poles = (),
                zeros = (),
                poles_overlay = fitter_factor.poles * fitter_factor.poles_overlay,
                zeros_overlay = fitter_factor.zeros * fitter_factor.zeros_overlay,
                delay_s = self.fitter.delay_s - fitter_factor.delay_s,
            )
            self.fitter = self.fitter.regenerate(ZPKrep = xrep)
        return

    def factorization_pop(self):
        factorization = self._fitter_factors.pop()
        xrep = representations.ZPKwData(
            ZPKrep = factorization.ZPKrep,
            zeros = self.fitter.zeros * factorization.zeros,
            poles = self.fitter.poles * factorization.poles,
            gain  = self.fitter.gain  * factorization.gain,
        )
        self.fitter = self.fitter.regenerate(ZPKrep = xrep)
        #TODO, check that the representation is preserved

    @contextlib.contextmanager
    def factorization(self, *args, **kwargs):
        self.factorization_push(*args, **kwargs)
        yield
        self.factorization_pop()

    def checkpoint(self, checkpoint):
        metadata = wavestate.bunch.Bunch()
        metadata.name = checkpoint
        metadata.log_idx           = self.log_number
        metadata.fitter_idx        = self(self._fitters) - 1
        metadata.fitter_current    = self._fitter_current
        metadata.fitter_lowres_avg = self._fitter_lowres_avg
        metadata.fitter_lowres_max = self._fitter_lowres_max
        metadata.fitter_lowres_med = self._fitter_lowres_med
        metadata.fitter            = self.fitter.copy()
        metadata.fitter_special    = wavestate.bunch.Bunch()
        for name, fitter in self.fitter_special.items():
            metadata.fitter_special[name] = fitter.copy()
        metadata.fitter_factors = list(self._fitter_factors)
        self._checkpoints.append(metadata)

    def invalidate_fitters(self):
        """
        If an algorithm a
        TODO
        """
        for fmeta in self._fitters:
            fmeta.valid = False
        return

    def log(
            self,
            *args,
            **kwargs
    ):
        """
        First argument is the level, should include a log group, which must be one of
        ['info', 'debug', 'warn', 'alert', 'rationale', 'progress']
        """
        level = args[0]
        if isinstance(level, int):
            level = args[0]
            args = args[1:]
            group = kwargs.setdefault('group', 'info')
        else:
            level = -1
            group = kwargs.setdefault('group', 'debug')
            #TODO print line and file upon hint request
            #args = args

        header = self.log_header_stack
        kwargs['header'] = header
        kwargs['time'] = time.time()
        kwargs['time_start'] = self.mtime_start

        if self.hint('log_off', default = False):
            return

        kwargs['args'] = args
        #TODO, merge if consecutive with the same parameters
        self._logs.append(
            kwargs
        )
        self.log_number += 1

        #FOR LIVE PRINTING
        if group == 'info':
            log_mod_level = logging.INFO
            group_character = 'I'
            level_limit = self.hint([
                'log_level_info',
                'log_level',
            ], default = 8)
        elif group == 'debug':
            log_mod_level = logging.DEBUG
            group_character = 'D'
            level_limit = self.hint([
                'log_level_debug',
                'log_level',
            ], default = 8)
        elif group == 'warn':
            log_mod_level = logging.WARNING
            group_character = 'W'
            level_limit = self.hint([
                'log_level_warn',
                'log_level',
            ], default = 8)
        elif group == 'alert':
            log_mod_level = logging.WARNING
            group_character = 'A'
            level_limit = self.hint([
                'log_level_alert',
                'log_level',
            ], default = 8)
        elif group == 'rationale':
            log_mod_level = logging.INFO
            group_character = 'R'
            level_limit = self.hint([
                'log_level_rationale',
                'log_level',
            ], default = 8)
        elif group == 'progress':
            log_mod_level = logging.INFO
            group_character = 'P'
            level_limit = self.hint([
                'log_level_progress',
                'log_level',
            ], default = 8)
        else:
            raise RuntimeError("Unrecognized log grouping")

        if self.hint('log_print', default = True) and level <= level_limit:
            hint_log_stdout = self.hint('log_stdout', default = True)
            if hint_log_stdout not in [None, True, False]:
                lfile = hint_log_stdout
            else:
                lfile = sys.stdout

            header = self.log_header_stack
            header_len = len(header)

            prefix = "{}{} {: >6.2f} {}".format(
                level if level >= 0 else '-',
                group_character,
                kwargs['time'] - kwargs['time_start'],
                '  ' * header_len
            )

            #TODO, make these take a header argument
            if not self.hint('logging_use', default = False):
                def pfunc(*args, **kwargs):
                    print(*args, **kwargs)
            else:
                def pfunc(*args, **kwargs):
                    kwargs.pop('file', None)
                    logging.log(log_mod_level + 9 - level, *args, **kwargs)

            if header_len > self.log_header_printed:
                pfunc(
                    "{}:{}:".format(
                        '-' * (len(prefix)),
                        ":".join(header)
                    ),
                    file = lfile
                )
                self.log_header_printed = header_len
                #tag that the header has been printed

            hint_log_stderr = self.hint('log_stderr', default = True)
            if hint_log_stderr and group == 'warn':
                if hint_log_stderr not in [None, True, False]:
                    lfile = hint_log_stderr
                else:
                    lfile = sys.stderr
            else:
                lfile = sys.stdout

            arg_lines = [[]]
            for arg in args:
                if isinstance(arg, (str, unicode)):
                    if '\n' in arg:
                        arg = padding_remove(arg)
                    arg_spl = arg.split('\n')
                    arg_lines[-1].append(arg_spl[0])
                    for subline in arg_spl[1:]:
                        arg_lines.append([subline])
                else:
                    arg_lines[-1].append(arg)

            #TODO, have pfunc do this splitting
            pfunc(
                prefix, *arg_lines[0],
                file = lfile
            )
            for argsl in arg_lines[1:]:
                pfunc(
                    ' ' * len(prefix), *argsl,
                    file = lfile
                )
        return

    def log_debug(self, *args, **kwargs):
        kwargs['group'] = 'debug'
        self.log(*args, **kwargs)

    def log_warn(self, *args, **kwargs):
        kwargs['group'] = 'warn'
        self.log(*args, **kwargs)

    def log_alert(self, *args, **kwargs):
        kwargs['group'] = 'alert'
        self.log(*args, **kwargs)

    def log_info(self, *args, **kwargs):
        kwargs['group'] = 'info'
        self.log(*args, **kwargs)

    def log_rationale(self, *args, **kwargs):
        kwargs['group'] = 'rationale'
        self.log(*args, **kwargs)

    def log_progress(self, *args, **kwargs):
        kwargs['group'] = 'progress'
        self.log(*args, **kwargs)

    @contextlib.contextmanager
    def log_heading(self, header):
        save_stack = self.log_header_stack
        self.log_header_stack = save_stack + (header,)
        #TODO, auto print header on command?
        yield
        self.log_header_stack = save_stack
        if self.log_header_printed > len(save_stack):
            self.log_header_printed = len(save_stack)
