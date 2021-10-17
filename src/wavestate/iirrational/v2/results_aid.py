#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""

import contextlib
import numpy as np

from .arguments import ArgumentError
from .. import fitters_ZPK
import warnings
import collections

from ..TFmath import match_SOS_pairs
from .. import representations

class ResultsAid(object):

    def __init__(self, fit_aid, kw = None):
        self._current_order = None
        self._fitter        = None
        self.kw             = kw

        self.initial_arguments = {}
        self.fitting_function = None

        self.fit_aid = fit_aid

        self.__fitters_by_order = None
        self._extra_data_by_order = dict()
        return

    @property
    def fitter(self):
        return self._fitter

    @property
    def _fitter_extra(self):
        extra = self._extra_data_by_order.setdefault(self._current_order, {})
        return extra

    @property
    def fitterRI(self):
        extra = self._fitter_extra
        fitterRI = extra.get('fitterRI', None)
        if fitterRI is None:
            fitterRI = self.fitter.regenerate(
                coding_map = fitters_ZPK.coding_maps.RI
            )
            extra['fitterRI'] = fitterRI
        return fitterRI

    @property
    def _fitters_by_order(self):
        """
        Returns four arrays, an array of fitter orders, residuals, the
        fitters themselves, and a bool array indicating if each fitter is the
        local best (the fitter of lowest residuals with equal or lower order).
        They are sorted by the fitter order such that binary
        search may be performed on the order array.
        """
        if self.__fitters_by_order is None:
            fitters_by_order = {}
            for fmeta in self.fit_aid._fitters:
                if not fmeta.valid:
                    continue

                fgroup = fitters_by_order.setdefault(fmeta.fitter.order, [])
                fgroup.append(fmeta.fitter)
            sorted_orders = []
            sorted_fitters = []
            sorted_residuals = []
            sorted_local_best = []
            lowest_res = float('inf')
            for order, fgroup in sorted(fitters_by_order.items()):
                fgroup.sort(key = lambda f: f.residuals_average)
                sorted_orders.append(order)
                sorted_fitters.append(fgroup[0])
                resavg = fgroup[0].residuals_average
                sorted_residuals.append(resavg)
                if resavg < lowest_res:
                    lowest_res = resavg
                    sorted_local_best.append(True)
                else:
                    sorted_local_best.append(False)
            sorted_local_best =np.asarray(sorted_local_best)
            self.__fitters_by_order = (
                np.asarray(sorted_orders),
                np.asarray(sorted_residuals),
                # currently a bug in declarative that attribute misses cause RuntimeError rather than AttributeError
                #numpy now duck-type checks this, and the wrong exception is failing. Instead, casting to list
                list(sorted_fitters),
                sorted_local_best,
            )
        return self.__fitters_by_order

    @property
    def order(self):
        return self._current_order

    def residuals_by_order(self):
        obo, rbo, fbo, lbo  = self._fitters_by_order
        return obo, rbo

    def choose(self, order = None, direction = None, **kwargs):
        """
        Chooses the filter with the lowest residual of order equal or less than
        requested. The data2filter algorithm will already have chosen a filter,
        but this allows investigating alternatives among the ones explored by
        the algorithm.

        if order is None, then the order is taken to be the current order.

        Can take argument "direction" to specify which nearest order to find
         - "<=" or "le": less than or equal to the specified order
         - ">=" or "ge": greater than or equal to the specified order
         - "<" or "l": less than the specified order
         - ">" or "g": greater than or equal to the specified order

        Note that calling
        results.choose(direction = '>') will return the fitter with the next
        order above the current one, and calling results.choose(direction = '<')
        gives the next lower order. So that one can investigate the nearest
        alternative choices.
        """
        if direction is None:
            direction = '<='

        if order is None:
            order = self.order

        dmap = {
            'ge' : '>=',
            'g' :  '>',
            '>=' : '>=',
            '>' :  '>',
            'le' : '<=',
            'l' :  '<',
            '<=' : '<=',
            '<' :  '<',
        }

        try:
            direction = dmap[direction]
        except KeyError:
            raise ArgumentError((
                "'direction' argument must be one of {}"
            ).format(dmap.keys))

        #using the while just as a goto statement with 'break' keyword
        while True:
            obo, rbo, fbo, lbo = self._fitters_by_order
            idx_ord = np.searchsorted(obo, order)
            if idx_ord == 0:
                if direction in ('<=', '<'):
                    break
            elif idx_ord >= len(obo):
                if direction in ('>=', '>'):
                    idx_ord = len(obo) - 1
                    break
                idx_ord = len(obo) - 1

            if order == obo[idx_ord]:
                if direction == '<':
                    idx_ord -= 1
                elif direction == '>':
                    idx_ord += 1

            if direction in ('<', '<='):
                idx_ord = np.argmin(rbo[:idx_ord + 1])
            break
        fitter = fbo[idx_ord]

        self._fitter = fitter
        self._current_order = obo[idx_ord]
        return fitter

    @contextlib.contextmanager
    def choice(self, order = None, direction = None, **kwargs):
        """
        Allows a temporary choice of fitter, for investigations without having
        to store the optimal choice. Takes the same arguments as choose().
        """
        fitter_prev = self._fitter
        self.choose(order = order, direction = direction, **kwargs)
        yield
        self._fitter = fitter_prev
        return

    def as_scipy_signal_ZPKsw(self, with_delay = False):
        """
        Returns the chosen filter in scipy format
        """
        #TODO
        from ..LIGO.conversion import filter2sortedZPK
        z, p, k = filter2sortedZPK(self.fitter)
        z = np.asarray(z) * 2 * np.pi
        p = np.asarray(p) * 2 * np.pi
        k = k * (2 * np.pi)**(len(p) - len(z))

        if not with_delay:
            if self.fitter.delay_s != 0:
                warnings.warn((
                    "fit export not including delay,"
                    " but filter has delay {}s. Call results.as_scipy_signal_ZPKsw(with_delay = True)"
                    " and the output will be ZPKD"
                ).format(self.fitter.delay_s))
            return z, p, k
        else:
            return z, p, k, self.fitter.delay_s

    def as_fitterMRF_ZPKz(
        self,
        t_sample_s = None,
        F_sample_Hz = None,
        F_nyquist_Hz = None,
    ):
        """
        Returns the chosen filter in scipy format for the z domain
        (must specify t_sample_s, or F_sample_Hz or F_nyquist_Hz). This
        transforms into the Z domain, then reruns the optimization to ensure a
        best fit in the transformed coordinates.
        """
        F_nyquist_Hz = argF_nyquist_Hz(
            t_sample_s   = t_sample_s,
            F_sample_Hz  = F_sample_Hz,
            F_nyquist_Hz = F_nyquist_Hz,
        )
        fitter_z = self._fitter_extra.get(F_nyquist_Hz, None)
        if fitter_z is None:
            xrep = self.fitter.ZPKrep.change_F_Nyquist_Hz(F_nyquist_Hz)
            raise NotImplementedError()
            self._fitter_extra[F_nyquist_Hz] = fitter_z
        return fitter_z

    def as_scipy_signal_ZPKz(
        self,
        t_sample_s = None,
        F_sample_Hz = None,
        F_nyquist_Hz = None,
    ):
        """
        Returns the chosen filter in scipy format for the z domain
        (must specify t_sample_s, or F_sample_Hz or F_nyquist_Hz). This
        transforms into the Z domain, then reruns the optimization to ensure a
        best fit in the transformed coordinates.
        """
        fitter_z = self.as_fitterMRF_ZPKz(
            t_sample_s   = t_sample_s,
            F_sample_Hz  = F_sample_Hz,
            F_nyquist_Hz = F_nyquist_Hz,
        )
        raise NotImplementedError()
        return

    def as_ZPKrep(
        self,
    ):
        """
        Returns the filter as an IIRrational ZPKwData object, for use with other
        fitters in the IIRrational package.
        """
        return self.fitter.ZPKrep

    @property
    def ZPKrep(self):
        return self.as_ZPKrep()

    def as_ZPKTF(
        self, delay2ZPK_args
    ):
        """
        Returns the filter as an IIRrational ZPKTF object, for use
        with mathematical analysis. Note, if the delay is nonzero,
        some conversion of the delay into ZPK is required. If the delay is
        positive (as is usual), the some order of RHP-zero/LHP-pole pairs must
        be added to approximate the delay.
        """
        #TODO
        raise NotImplementedError()

    def as_matlab_str_ZPKsw(
        self,
    ):
        """
        Returns the chosen filter in matlab format for the z domain
        (must specify t_sample_s, or F_sample_Hz or F_nyquist_Hz)
        """
        #TODO
        raise NotImplementedError()

    def as_foton_ZPKsf(
        self,
        with_delay = False,
    ):
        """
        Returns the chosen filter as a string foton formatting. Uses the
        frequency domain representation
        """
        from ..LIGO.conversion import filter2fotonZPK
        ZPK = filter2fotonZPK(
            self.fitter.ZPKrep,
            plane = 'f',
            zpk_output = True,
        )
        if not with_delay:
            if self.fitter.delay_s != 0:
                warnings.warn((
                    "fit export not including delay,"
                    " but filter has delay {}s. Call results.as_foton_ZPKsf(with_delay = True)"
                    " and the output will be ZPKD"
                ).format(self.fitter.delay_s))
            return ZPK
        else:
            return tuple(ZPK) + (self.fitter.delay_s,)

    def as_foton_ZPKsn(
        self,
        with_delay = False,
    ):
        """
        Returns the chosen filter as a string foton formatting. Uses the
        frequency domain representation
        """
        from ..LIGO.conversion import filter2fotonZPK
        ZPK = filter2fotonZPK(
            self.fitter.ZPKrep,
            plane = 'n',
            zpk_output = True,
        )
        if not with_delay:
            if self.fitter.delay_s != 0:
                warnings.warn((
                    "fit export not including delay,"
                    " but filter has delay {}s. Call results.as_foton_ZPKsn(with_delay = True)"
                    " and the output will be ZPKD"
                ).format(self.fitter.delay_s))
            return ZPK
        else:
            return tuple(ZPK) + (self.fitter.delay_s,)

    def as_refine_str(self, sigfigs = 5):
        zzpp_pairs = match_SOS_pairs(
            self.fitter.zeros.r, self.fitter.zeros.c,
            self.fitter.poles.r, self.fitter.poles.c,
            F_nyquist_Hz = None
        )
        zs = []
        ps = []

        def append(rs, r):
            if r is not None:
                if r.imag != 0:
                    if r.imag > 0:
                        rs.append("{real:.{sf}f}+{imag:.{sf}f}j".format(real = r.real, imag = r.imag, sf = sigfigs))
                    else:
                        rs.append("{real:.{sf}f}-{imag:.{sf}f}j".format(real = r.real, imag = -r.imag, sf = sigfigs))
                else:
                    rs.append("{real:.{sf}f}".format(real = r.real, sf = sigfigs))

        for (z1, z2, p1, p2) in zzpp_pairs:
            append(zs, z1)
            append(zs, z2)
            append(ps, p1)
            append(ps, p2)
        return "--zeros={} --poles={}".format(",".join(zs), ",".join(ps))

    def as_foton_str_ZPKsf(
        self,
        annotate_pairs = False,
    ):
        """
        Returns the chosen filter as a string foton formatting. Uses the
        frequency domain representation
        """
        from ..LIGO.conversion import filter2fotonZPK
        return filter2fotonZPK(
            self.fitter.ZPKrep,
            annotate_pairs = annotate_pairs,
            plane = 'f',
        )

    def as_foton_str_ZPKsn(
        self,
        annotate_pairs = False,
    ):
        """
        Returns the chosen filter as a string foton formating. Uses the
        "normalized" frequency domain representation.
        """
        #TODO
        from ..LIGO.conversion import filter2fotonZPK
        return filter2fotonZPK(
            self.fitter.ZPKrep,
            annotate_pairs = annotate_pairs,
            plane = 'n',
        )

    def as_foton_str_ZPKz(
        self,
        t_sample_s = None,
        F_sample_Hz = None,
        F_nyquist_Hz = None,
    ):
        """
        Returns the chosen filter in foton format native for the z domain
        (must specify t_sample_s, or F_sample_Hz or F_nyquist_Hz). This
        transforms into the Z domain, then reruns the optimization to ensure a
        best fit in the transformed coordinates.
        """
        #TODO
        from ..LIGO.conversion import filter2fotonZPK
        fitter_z = self.as_fitterMRF_ZPKz(
            t_sample_s   = t_sample_s,
            F_sample_Hz  = F_sample_Hz,
            F_nyquist_Hz = F_nyquist_Hz,
        )
        return filter2fotonZPK(fitter_z.ZPKrep, plane = 'z')

    def _fitter_export_dict(self, fitter):
        from ..LIGO.conversion import filter2fotonZPK
        zpk = tuple(fitter.ZPKrep.ZPK)
        d = collections.OrderedDict([
            ("zpk", collections.OrderedDict([
                ("z" , zpk[0]),
                ("p" , zpk[1]),
                ("k" , zpk[2]),
            ])),
            ("delay_s", fitter.ZPKrep.delay_s),
            ("gain", fitter.gain),
            ("poles", collections.OrderedDict([
                ("r" , fitter.poles.r),
                ("c" , fitter.poles.c),
            ])),
            ("zeros", collections.OrderedDict([
                ("r" , fitter.zeros.r),
                ("c" , fitter.zeros.c),
            ])),
            ("poles_overlay", collections.OrderedDict([
                ("r" , fitter.poles_overlay.r),
                ("c" , fitter.poles_overlay.c),
            ])),
            ("zeros_overlay", collections.OrderedDict([
                ("r" , fitter.zeros_overlay.r),
                ("c" , fitter.zeros_overlay.c),
            ])),
            ("residuals", fitter.residuals_average),
            ("residuals_max", fitter.residuals_max),
            ("residuals_med", fitter.residuals_med),
            ("foton_Sf", filter2fotonZPK(
                fitter.ZPKrep,
                annotate_pairs = False,
                plane = 'f',
            ))
        ])
        return d

    def choice_export_dict(self):
        return self._fitter_export_dict(self.fitter)

    def alternatives_export_dict(self, local_best = True):
        obo, rbo, fbo, lbo  = self._fitters_by_order

        #only returns or displays those which are the best for that order or below
        if local_best:
            obo = obo[lbo]
            rbo = rbo[lbo]
            fbo = [fbo[idx] for idx, val in enumerate(lbo) if val]

        exd = collections.OrderedDict()
        for order, fit in zip(obo, fbo):
            key = 'order_{:0>2d}'.format(order)
            exd[key] = self._fitter_export_dict(fit)

        return exd


def argF_nyquist_Hz(
    t_sample_s = None,
    F_sample_Hz = None,
    F_nyquist_Hz = None,
):
    if t_sample_s is not None:
        _F_nyquist_Hz = 1/(2 * t_sample_s)
        if F_nyquist_Hz is not None and _F_nyquist_Hz != F_nyquist_Hz:
            raise ArgumentError(
                "t_sample_s and F_nyquist_Hz both"
                " given, but not consistent"
            )
        F_nyquist_Hz = _F_nyquist_Hz
    if F_sample_Hz is not None:
        _F_nyquist_Hz = F_sample_Hz/2
        if F_nyquist_Hz is not None and _F_nyquist_Hz != F_nyquist_Hz:
            raise ArgumentError(
                "F_sample_Hz and F_nyquist_Hz both"
                " given, but not consistent"
            )
        F_nyquist_Hz = _F_nyquist_Hz
    return F_nyquist_Hz

