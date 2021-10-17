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

from declarative import (
    depB_property,
    DepBunch,
)
import numpy as np
from ... import TFmath


def abs_sq(x):
    return x.real**2 + x.imag**2


#root_constraints:
class RootConstraints(object):
    no_constraint = frozenset()  # real Even real Odd
    mirror_real = frozenset(['mirror_real'])  # real Even real Odd
    mirror_imag = frozenset(['mirror_imag'])  # real Even imaginary Odd
    mirror_disc = frozenset(['mirror_disc'])  # palendromic standard polynomials

    mirror_quad = mirror_real | mirror_imag

    mirror_real_disc = mirror_real | mirror_disc
    mirror_imag_disc = mirror_imag | mirror_disc
    mirror_quad_disc = mirror_quad | mirror_disc


root_constraints = RootConstraints()

#rootBunches must have the following elements
#rB.constraint which defines how to interpret
#rB.z is the number of zero roots
#rB.c are complex roots
#rB.r are real roots
#rB.i are imaginary roots
#rB.dc are disc complex roots
#rB.dr are disc real roots
#rB.di are disc imaginary roots
#rB.u are unsorted roots (not allowed for many applications)

#if the mirror_imag is not specified, imaginary roots will not be separated
#if the mirror_real is not specified, real roots will not be separated
#if the mirror_disc is not specified, disc roots will not be separated

#if mirror flags are given, only roots non-redundant are provided in the lists
#zero roots will always be in the real list
#

class RootBunch(DepBunch):
    root_constraints = root_constraints
    def __init__(
        self,
        constraint = None,
        z  = None,
        c  = None,
        r  = None,
        i  = None,
        dc = None,
        dr = None,
        di = None,
        u  = (),
        copy = None,
        **kwargs
    ):
        if copy is not None:
            #TODO, make this more concise
            return super(RootBunch, self).__init__(
                self,
                constraint = constraint,
                z  = z, c  = i,
                r  = r, i  = i,
                dc = dc, dr = dr,
                di = di, u  = u,
                copy = copy,
                **kwargs
            )
        remainder = constraint - root_constraints.mirror_quad_disc
        if remainder:
            raise RuntimeError("Unknown root_constraints: {0}".format(remainder))
        if copy is not None:
            super(RootBunch, self).__init__(copy = copy, **kwargs)
            return

        if constraint is None:
            raise TypeError("Must specify constraint argument")

        super(RootBunch, self).__init__(**kwargs)
        self.constraint = constraint
        if z is not None:
            self.z  = int(z)
        if c is not None:
            self.c  = np.asarray(c)
        if r is not None:
            #TODO, check realness
            self.r  = np.asarray(r)
        if i is not None:
            #TODO, check imaginaryness
            self.i  = np.asarray(i)
        if dc is not None:
            self.dc = np.asarray(dc)
        if dr is not None:
            #TODO, check realness
            self.dr = np.asarray(dr)
        if di is not None:
            #TODO, check imaginaryness
            self.di = np.asarray(di)
        self.u  = np.asarray(u)

        if root_constraints.mirror_disc <= constraint:
            assert(dc is not None)
            if root_constraints.mirror_quad <= constraint:
                assert(z is not None)
                assert(c is not None)
                assert(r is not None)
                assert(i is not None)
                assert(dc is not None)
                assert(dr is not None)
                assert(di is not None)
                #TODO, needs even more asserts
                assert(np.all(self.r.real < 0))
                assert(np.all(self.i.real > 0))
                assert(np.all(self.r.imag == 0))
                assert(np.all(self.i.imag == 0))
                assert(np.all(self.dr.real < 0))
                assert(np.all(self.di.real > 0))
                assert(np.all(self.dr.imag == 0))
                assert(np.all(self.di.imag == 0))
            elif root_constraints.mirror_real <= constraint:
                assert(z is None)
                assert(r is not None)
                assert(c is not None)
                assert(dr is not None)
                assert(dc is not None)
                #TODO, needs even more asserts
                assert(np.all(self.c.imag > 0))
                assert(np.all(self.dc.imag > 0))
                assert(np.all(self.r.imag == 0))
                assert(np.all(self.dr.imag == 0))
            elif root_constraints.mirror_imag <= constraint:
                assert(z is None)
                assert(i is not None)
                assert(c is not None)
                assert(di is not None)
                assert(dc is not None)
                #TODO, needs even more asserts
                assert(np.all(self.c.real > 0))
                assert(np.all(self.dc.real > 0))
                assert(np.all(self.i.imag == 0))
                assert(np.all(self.di.imag == 0))
        else:
            assert(dc is None)
            assert(dr is None)
            assert(di is None)
            if root_constraints.mirror_quad <= constraint:
                assert(z is not None)
                assert(c is not None)
                assert(r is not None)
                assert(i is not None)
                assert(np.all(self.c.real < 0))
                assert(np.all(self.c.imag > 0))
                assert(np.all(self.r.real < 0))
                assert(np.all(self.i.real > 0))
                assert(np.all(self.r.imag == 0))
                assert(np.all(self.i.imag == 0))
            elif root_constraints.mirror_real <= constraint:
                assert(z is None)
                assert(r is not None)
                assert(c is not None)
                assert(i is None)
                assert(np.all(self.c.imag > 0))
                assert(np.all(self.r.imag == 0))
            elif root_constraints.mirror_imag <= constraint:
                assert(z is None)
                assert(i is not None)
                assert(c is not None)
                assert(r is None)
                assert(np.all(self.c.real > 0))
                assert(np.all(self.i.imag == 0))
            else:
                assert(z is None)
                assert(i is None)
                assert(c is None)
                assert(r is None)
        return

    def __len__(self):
        return self.length

    def clear(self):
        self._count = self._count + 1

    def root_str_list(self):
        slist = []
        if root_constraints.mirror_disc <= self.constraint:
            raise NotImplementedError()
            if root_constraints.mirror_quad <= self.constraint:
                pass
            elif root_constraints.mirror_real <= self.constraint:
                pass
            elif root_constraints.mirror_imag <= self.constraint:
                pass
            else:
                #palidromic but complex
                pass
        else:
            if root_constraints.mirror_quad <= self.constraint:
                for r in sorted(self.r):
                    slist.append("±{}".format(r))
                for i in sorted(self.i):
                    slist.append("±{}j".format(i))
                for c in sorted(self.c, key = lambda v : v.imag):
                    slist.append("±{}±{}j".format(c.real, c.imag))
            elif root_constraints.mirror_real <= self.constraint:
                for r in sorted(self.r):
                    slist.append("{}".format(r))
                for c in sorted(self.c, key = lambda v : v.imag):
                    slist.append("{}±{}j".format(c.real, c.imag))
            elif root_constraints.mirror_imag <= self.constraint:
                for i in sorted(self.i):
                    slist.append("{}j".format(i))
                for c in sorted(self.c, key = lambda v : v.imag):
                    if c.imag == 0:
                        slist.append("±{}".format(c.real))
                    elif c.imag > 0:
                        slist.append("±{}+{}j".format(c.real, c.imag))
                    else:
                        slist.append("±{}-{}j".format(c.real, -c.imag))
            else:
                pass
        for u in sorted(self.u, key = lambda v : v.real):
            if u.imag == 0:
                slist.append("{}".format(u.real))
            elif u.imag > 0:
                slist.append("{}+{}j".format(u.real, u.imag))
            else:
                slist.append("{}-{}j".format(u.real, -u.imag))
        return slist

    def __str__(self):
        slist = self.root_str_list()
        #figure out how to make this pretty print with lines!

        s = "rB({})".format(', '.join(slist))
        #TODO, make this nicer
        if sys.version_info < (3, 0):
            s = s.encode('utf-8')
        return s

    # TODO make a REPR

    @depB_property
    def _count(self, val = 0):
        return val

    @depB_property
    def length(self):
        self.dependencies('length')
        if root_constraints.mirror_disc <= self.constraint:
            raise NotImplementedError()
            if root_constraints.mirror_quad <= self.constraint:
                pass
            elif root_constraints.mirror_real <= self.constraint:
                pass
            elif root_constraints.mirror_imag <= self.constraint:
                pass
            else:
                #palidromic but complex
                pass
        else:
            if root_constraints.mirror_quad <= self.constraint:
                return (self.z
                        + 2 * self.r.size
                        + 2 * self.i.size
                        + 4 * self.c.size
                        + self.u.size)
            elif root_constraints.mirror_real <= self.constraint:
                return(
                    self.r.size
                    + 2 * self.c.size
                    + self.u.size
                )
            elif root_constraints.mirror_imag <= self.constraint:
                return (
                    self.i.size
                    + 2 * self.c.size
                    + self.u.size
                )
            else:
                return len(self.u)

    @depB_property
    def fullplane(self):
        self.dependencies('length')
        if root_constraints.mirror_disc <= self.constraint:
            raise NotImplementedError()
            if root_constraints.mirror_quad <= self.constraint:
                pass
            elif root_constraints.mirror_real <= self.constraint:
                pass
            elif root_constraints.mirror_imag <= self.constraint:
                pass
            else:
                #palidromic but complex
                pass
        else:
            if root_constraints.mirror_quad <= self.constraint:
                a = np.empty((
                    self.z
                    + 2 * self.r.size
                    + 2 * self.i.size
                    + 4 * self.c.size
                    + self.u.size,
                ), dtype=self.c.dtype)
                a[:self.z] = 0
                offset = self.z
                offset_2 = offset + 2 * self.r.size
                v = a[offset:offset_2]
                v[::2] = self.r
                v[1::2] = -self.r
                offset = offset_2
                offset_2 = offset + 2 * self.i.size
                v = a[offset:offset_2]
                v[::2] = 1j*self.i
                v[1::2] = -1j*self.i
                offset = offset_2
                offset_2 = offset + 4 * self.c.size
                v = a[offset:offset_2]
                v[::4] = self.c
                v[1::4] = self.c.conjugate()
                v[2::4] = -self.c
                v[3::4] = -self.c.conjugate()
                return a
                offset = offset_2
                offset_2 = offset + self.u.size
                v = a[offset:offset_2]
                v[:] = self.u
            elif root_constraints.mirror_real <= self.constraint:
                a = np.empty((
                    self.r.size
                    + 2 * self.c.size
                    + self.u.size
                ), dtype=self.c.dtype)
                offset = 0
                offset_2 = offset + 1 * self.r.size
                v = a[offset:offset_2]
                v[:] = self.r
                offset = offset_2
                offset_2 = offset + 2 * self.c.size
                v = a[offset:offset_2]
                v[::2] = self.c
                v[1::2] = self.c.conjugate()
                offset = offset_2
                offset_2 = offset + self.u.size
                v = a[offset:offset_2]
                v[:] = self.u
                return a
            elif root_constraints.mirror_imag <= self.constraint:
                a = np.empty((
                    self.i.size
                    + 2 * self.c.size
                    + self.u.size
                ), dtype=self.c.dtype)
                offset = 0
                offset_2 = offset + self.i.size
                v = a[offset:offset_2]
                v[:] = 1j*self.i
                offset = offset_2
                offset_2 = offset + 2 * self.c.size
                v = a[offset:offset_2]
                v[::2] = self.c
                v[1::2] = -self.c
                offset = offset_2
                offset_2 = offset + self.u.size
                v = a[offset:offset_2]
                v[:] = self.u
                return a
            else:
                return self.u

    def multiply_by(self, scale_factor):
        if np.imag(scale_factor) == 0:
            if root_constraints.mirror_disc <= self.constraint:
                raise NotImplementedError()
            else:
                if root_constraints.mirror_quad == self.constraint:
                    return self.__class__(
                        constraint = self.constraint,
                        z = self.z,
                        u = scale_factor * self.u,
                        c = scale_factor * self.c,
                        r = scale_factor * self.r,
                        i = scale_factor * self.i,
                    )
                elif root_constraints.mirror_real == self.constraint:
                    return self.__class__(
                        constraint = self.constraint,
                        u = scale_factor * self.u,
                        c = scale_factor * self.c,
                        r = scale_factor * self.r,
                    )
                elif root_constraints.mirror_imag == self.constraint:
                    return self.__class__(
                        constraint = self.constraint,
                        u = scale_factor * self.u,
                        c = scale_factor * self.c,
                        i = scale_factor * self.i,
                    )
                elif root_constraints.no_constraint == self.constraint:
                    return self.__class__(
                        constraint = self.constraint,
                        u = scale_factor * self.u,
                    )
                else:
                    raise RuntimeError("Unrecognized Root_constraints")

        elif np.real(scale_factor) == 0:
            scale_factor_i = scale_factor.imag
            if root_constraints.mirror_disc <= self.constraint:
                raise NotImplementedError()
            else:
                #TODO
                if root_constraints.mirror_quad == self.constraint:
                    scale_factor_abs = abs(scale_factor_i)
                    return self.__class__(
                        constraint = root_constraints.mirror_quad,
                        z = self.z,
                        u = 1j * scale_factor_i * self.u,
                        c = (1j * scale_factor_abs * self.c).conjugate(),
                        r = -scale_factor_abs * self.i,
                        i = -scale_factor_abs * self.r,
                    )
                elif root_constraints.mirror_real == self.constraint:
                    if scale_factor_i > 0:
                        return self.__class__(
                            constraint = root_constraints.mirror_imag,
                            u = scale_factor_i * self.u,
                            c = (1j * scale_factor_i * self.c.conjugate()),
                            i = scale_factor_i * self.r,
                        )
                    else:
                        return self.__class__(
                            constraint = root_constraints.mirror_imag,
                            u = scale_factor_i * self.u,
                            c = (1j * scale_factor_i * self.c),
                            i = scale_factor_i * self.r,
                        )
                elif root_constraints.mirror_imag == self.constraint:
                    if scale_factor_i > 0:
                        return self.__class__(
                            constraint = root_constraints.mirror_real,
                            u = scale_factor_i * self.u,
                            c = (1j * scale_factor_i * self.c),
                            r = -scale_factor_i * self.i,
                        )
                    else:
                        return self.__class__(
                            constraint = root_constraints.mirror_real,
                            u = scale_factor_i * self.u,
                            c = (1j * scale_factor_i * self.c.conjugate()),
                            r = -scale_factor_i * self.i,
                        )
                elif root_constraints.no_constraint == self.constraint:
                    return self.__class__(
                        constraint = self.constraint,
                        u = scale_factor_i * self.u,
                    )
                else:
                    raise RuntimeError("Unrecognized Root_constraints")
        else:
            raise NotImplementedError("Cannot rotate continuous angles")

    def _multiply_RB(self, other):
        """

        """
        if other.constraint != self.constraint:
            raise RuntimeError("Can't multiply root bunches with differing constraints")

        mR = bool(root_constraints.mirror_real & self.constraint)
        mI = bool(root_constraints.mirror_imag & self.constraint)
        mD = bool(root_constraints.mirror_disc & self.constraint)
        kw = dict()
        kw['u'] = np.concatenate([self.u, other.u])
        if mR or mI:
            kw['c'] = np.concatenate([self.c, other.c])
            if mR:
                kw['r'] = np.concatenate([self.r, other.r])
                if mI:
                    kw['i'] = np.concatenate([self.i, other.i])
                    kw['z'] = self.z + other.z
                    if mD:
                        kw['di'] = np.concatenate([self.di, other.di])
                        kw['dr'] = np.concatenate([self.dr, other.dr])
                        kw['dc'] = np.concatenate([self.dc, other.dc])
                elif mD:
                    kw['dr'] = np.concatenate([self.dr, other.dr])
                    kw['dc'] = np.concatenate([self.dc, other.dc])
            else:
                #mR not true, mI true
                kw['i'] = np.concatenate([self.i, other.i])
                if mD:
                    kw['di'] = np.concatenate([self.di, other.di])
                    kw['dc'] = np.concatenate([self.dc, other.dc])
        return self.__class__(
            constraint = self.constraint,
            **kw
        )

    def __mul__(self, other):
        if isinstance(other, RootBunch):
            return self._multiply_RB(other)
        return self.multiply_by(other)

    def __rmul__(self, other):
        return self.multiply_by(other)

    def __truediv__(self, other):
        return self.multiply_by(1/other)

    def __div__(self, other):
        return self.multiply_by(1/other)

    def val_lnG(self, X, h = 1, lnG = 0):
        """
        returns the value as if it were generated from a polynomial with last coefficient 1 given a coefficient representation
        and the X_scale. It computes this without having to convert to the coefficient representation (how efficient!).
        """
        X = np.asarray(X)
        h = np.array(h, copy = True, dtype = np.complex128)
        X, h = np.broadcast_arrays(X, h)

        #note that this modifies in-place
        def VfR(roots, h, lnG):
            if len(roots) == 0:
                return h, lnG
            roots = np.asarray(roots)
            mlen = len(roots)
            group_len = 5
            for idx in range((mlen - 1) // group_len + 1):
                r = roots[idx * group_len : (idx + 1) * group_len]
                h = h * np.polynomial.polynomial.polyvalfromroots(X, r)
                abs_max = np.max(abs_sq(h))**.5
                h /= abs_max
                lnG += np.log(abs_max)
            return h, lnG

        if len(self.u) > 0:
            h, lnG = VfR(self.u, h, lnG)

        mR = bool(root_constraints.mirror_real & self.constraint)
        mI = bool(root_constraints.mirror_imag & self.constraint)
        mD = bool(root_constraints.mirror_disc & self.constraint)
        if mR or mI:
            h, lnG = VfR(self.c, h, lnG)
            if mR:
                h, lnG = VfR(self.c.conjugate(), h, lnG)
                h, lnG = VfR(self.r, h, lnG)
                if mI:
                    h, lnG = VfR(-self.c, h, lnG)
                    h, lnG = VfR(-self.c.conjugate(), h, lnG)
                    h, lnG = VfR(-self.r, h, lnG)
                    h, lnG = VfR(1j*self.i, h, lnG)
                    h, lnG = VfR(-1j*self.i, h, lnG)
                    if mD:
                        h, lnG = VfR(1/self.c.conjugate(), h, lnG)
                        h, lnG = VfR(1/self.c, h, lnG)
                        h, lnG = VfR(-1/self.c.conjugate(), h, lnG)
                        h, lnG = VfR(-1/self.c, h, lnG)
                        h, lnG = VfR(1/self.r, h, lnG)
                        h, lnG = VfR(-1/self.r, h, lnG)
                        h, lnG = VfR(1j/self.i, h, lnG)
                        h, lnG = VfR(-1j/self.i, h, lnG)
                elif mD:
                    h, lnG = VfR(1/self.c.conjugate(), h, lnG)
                    h, lnG = VfR(1/self.c, h, lnG)
                    h, lnG = VfR(1/self.r, h, lnG)
            else:
                #mR not true, mI true
                h, lnG = VfR(-self.c, h, lnG)
                h, lnG = VfR(1j*self.i, h, lnG)
                if mD:
                    h, lnG = VfR(1/self.c.conjugate(), h, lnG)
                    h, lnG = VfR(-1/self.c.conjugate(), h, lnG)
                    h, lnG = VfR(1j/self.i, h, lnG)
        return h, lnG


class RBAlgorithms(object):
    root_constraints = root_constraints
    def __init__(
        self,
        strict    = True,
        line_tol  = 1e-6,
        match_tol = 1e-4,
        zero_tol  = 1e-10,
        lax_line_tol = 0,
        constraint_standard = None
    ):
        self.strict = strict
        self.lax_line_tol = lax_line_tol
        self.constraint_standard = constraint_standard

        def are_same(r1, r2):
            if abs(r1) < .8:
                return abs(r1 - r2) < match_tol
            else:
                return abs((r1/r2) - 1) < match_tol

        def are_real(r1):
            if abs(r1.real) < .8:
                return abs(np.imag(r1)) < line_tol
            else:
                return abs(np.imag(r1) / np.real(r1)) < line_tol

        def are_imag(r1):
            if abs(r1.imag) < .8:
                return abs(np.real(r1)) < line_tol
            else:
                return abs(np.real(r1) / np.imag(r1)) < line_tol

        def are_zero(r1):
            return abs(np.imag(r1)) < zero_tol

        self.are_same = np.vectorize(are_same, otypes = [bool])
        self.are_real = np.vectorize(are_real, otypes = [bool])
        self.are_imag = np.vectorize(are_imag, otypes = [bool])
        self.are_zero = np.vectorize(are_zero, otypes = [bool])

    def expect_atleast(
        self,
        rB,
        constraint = None,
        allow_unknown = False,
    ):
        if constraint is None:
            constraint = self.constraint_standard
            #root_constraints.no_constraint
        if not isinstance(rB, RootBunch):
            rB = RootBunch(
                u = np.asarray(rB),
                constraint = root_constraints.no_constraint,
            )

        if constraint <= rB.constraint:
            return rB

        return self.expect(
            rB,
            constraint = constraint,
            allow_unknown = allow_unknown
        )

    def expect(
        self,
        rB,
        constraint = None,
        allow_unknown = False,
    ):
        if constraint is None:
            constraint = self.constraint_standard
        #real root_constraints.mirrors first
        #imaginary root_constraints.mirrors second
        #now do disc root_constraints.mirrors
        remainder = constraint - root_constraints.mirror_quad_disc
        if remainder:
            raise RuntimeError("Unknown root_constraints: {0}".format(remainder))

        if not isinstance(rB, RootBunch):
            rB = RootBunch(
                u = np.asarray(rB),
                constraint = root_constraints.no_constraint,
            )
            return self.expect(
                rB,
                constraint,
                allow_unknown = allow_unknown
            )

        elif rB.constraint == constraint:
            if constraint == root_constraints.no_constraint:
                return rB
            elif not allow_unknown and len(rB.u) > 0:
                raise RuntimeError("Roots contain unpaired roots: {}".format(rB.u))
            return rB
        if rB.constraint <= constraint:
            return self._expect_elevate(
                rB, constraint,
                allow_unknown = allow_unknown,
            )
        elif constraint <= rB.constraint:
            return self._expect_reduce(
                rB, constraint,
                allow_unknown = allow_unknown,
            )
        else:
            #reduce to common root_constraints
            rB = self._expect_reduce(
                rB, (constraint & rB.constraint),
                allow_unknown = allow_unknown,
            )
            #now elevate to the new constraint
            return self._expect_elevate(
                rB, constraint,
                allow_unknown = allow_unknown,
            )
        raise NotImplementedError("Can't Get here")
        return

    def _expect_elevate(
        self,
        rB,
        constraint,
        allow_unknown
    ):
        #can only be called if rB.constraint <= constraint and rB.constraint != constraint
        assert(rB.constraint != constraint)
        def check_allow_unknown(rB):
            if not allow_unknown:
                if len(rB.u) > 0:
                    raise RuntimeError("Roots contain unpaired roots: {}".format(rB.u))

        def recurse(rB):
            check_allow_unknown(rB)
            if rB.constraint == constraint:
                return rB
            return self._expect_elevate(
                rB,
                constraint,
                allow_unknown = allow_unknown,
            )

        #elevate root_constraints
        if root_constraints.mirror_disc <= constraint:
            raise NotImplementedError("Can't generate root_constraints.mirror_disc root_constraints yet")
        elif root_constraints.mirror_real == rB.constraint:
            rB = self.MR2MQ(
                roots_c = rB.c,
                roots_r = rB.r,
                roots_u = rB.u,
            )
            return recurse(rB)
        elif root_constraints.mirror_imag == rB.constraint:
            rB = self.MI2MQ(
                roots_c = rB.c,
                roots_i = rB.i,
                roots_u = rB.u,
            )
            return recurse(rB)
        elif root_constraints.no_constraint == rB.constraint:
            if root_constraints.mirror_real <= constraint:
                rB = self.NC2MR(
                    roots_u = rB.u,
                )
                return recurse(rB)
            elif root_constraints.mirror_imag <= constraint:
                rB = self.NC2MI(
                    roots_u = rB.u,
                )
                return recurse(rB)
            else:
                raise NotImplementedError("Unknown final constraint")
        else:
            raise NotImplementedError("Unknown constraint in argument")
        return

    def _expect_reduce(
        self,
        rB,
        constraint,
        allow_unknown
    ):
        def check_allow_unknown(rB):
            if constraint == root_constraints.no_constraint:
                pass
            elif not allow_unknown and len(rB.u) > 0:
                    raise RuntimeError("Roots contain unpaired roots")
        def recurse(rB):
            check_allow_unknown(rB)
            return self._expect_reduce(
                rB,
                constraint,
                allow_unknown = allow_unknown,
            )

        if constraint == root_constraints.no_constraint:
            return RootBunch(
                u = rB.fullplane,
                constraint = root_constraints.no_constraint,
            )
        #reduce root_constraints
        raise NotImplementedError("Can not currently reduce root_constraints")
        if root_constraints.mirror_disc <= rB.constraint:
            raise NotImplementedError("Can't handle root_constraints.mirror disc root_constraints")
        elif root_constraints.mirror_quad <= rB.constraint:
            if root_constraints.mirror_real <= constraint:
                raise NotImplementedError()
                return
            else:
                #must be root_constraints.mirror_imag <= constraint
                raise NotImplementedError()
                return
            return
        elif root_constraints.mirror_real <= rB.constraint:
            raise NotImplementedError()
            return
        elif root_constraints.mirror_imag <= rB.constraint:
            raise NotImplementedError()
            return
        else:
            raise RuntimeError("Bad Logic")

    def NC2MR(self, roots_u):
        real_select = self.are_real(roots_u)
        roots_r = roots_u[real_select].real
        roots_u = roots_u[~real_select]

        pos_select = roots_u.imag > 0
        roots_c_neg = roots_u[~pos_select]
        roots_c_pos = roots_u[pos_select]
        rPB = TFmath.nearest_pairs(roots_c_pos, roots_c_neg.conjugate())
        if self.lax_line_tol > 0:
            roots_u = []
            roots_r2 = []
            def check_ins(u):
                if abs(u.real) > self.lax_line_tol:
                    if abs(u.imag) < self.lax_line_tol:
                        roots_r2.append(u.real)
                    else:
                        roots_u.append(u)
                else:
                    if abs(u.imag / u.real) < self.lax_line_tol:
                        roots_r2.append(u.real)
                    else:
                        roots_u.append(u)
            for u in rPB.l1_remain:
                check_ins(u)
            for u in rPB.l2_remain:
                check_ins(u.conjugate())
            roots_r = np.concatenate([roots_r, roots_r2])
        else:
            roots_u = list(rPB.l1_remain) + [r.conjugate() for r in rPB.l2_remain]
        roots_c = []
        for r1, r2 in rPB.r12_list:
            if not self.strict or self.are_same(r1, r2):
                #roots_c.append((r1 + r2) / 2)
                #TODO, this seems to work better, not clear why..
                roots_c.append(r1)
            else:
                roots_u.append(r1)
                roots_u.append(r2.conjugate())
        roots_c = np.array(roots_c)
        roots_u = np.array(roots_u)
        return RootBunch(
            constraint = root_constraints.mirror_real,
            c = roots_c,
            r = roots_r,
            u = roots_u,
        )

    def NC2MI(self, roots_u):
        imag_select = self.are_imag(roots_u)
        roots_i = roots_u[imag_select].imag
        roots_u = roots_u[~imag_select]

        pos_select = roots_u.real > 0
        roots_c_neg = roots_u[~pos_select]
        roots_c_pos = roots_u[pos_select]
        rPB = TFmath.nearest_pairs(roots_c_pos, -roots_c_neg.conjugate())

        #TODO, put this logic in the NC2MR

        #fill the list with as many pairs as possible
        r12_list_full = rPB.r12_list
        while rPB.r12_list:
            rPB = TFmath.nearest_pairs(rPB.l1_remain, rPB.l2_remain)
            r12_list_full.extend(rPB.r12_list)
        rPB.r12_list = r12_list_full

        if self.lax_line_tol > 0:
            roots_u = []
            roots_i2 = []
            def check_ins(u):
                if abs(u.imag) > self.lax_line_tol:
                    if abs(u.real) < self.lax_line_tol:
                        roots_i2.append(u.imag)
                    else:
                        roots_u.append(u)
                else:
                    if abs(u.real / u.imag) < self.lax_line_tol:
                        roots_i2.append(u.imag)
                    else:
                        roots_u.append(u)
            for u in rPB.l1_remain:
                check_ins(u)
            for u in rPB.l2_remain:
                check_ins(-u.conjugate())
            roots_i = np.concatenate([roots_i, roots_i2])
        else:
            roots_u = list(rPB.l1_remain) + [-r.conjugate() for r in rPB.l2_remain]

        roots_c = []
        for r1, r2 in rPB.r12_list:
            if not self.strict or self.are_same(r1, r2):
                #roots_c.append((r1 + r2) / 2)
                #TODO, this seems to work better, not clear why..
                roots_c.append(r1)
            else:
                roots_u.append(r1)
                roots_u.append(-r2.conjugate())
        roots_c = np.array(roots_c)
        roots_u = np.array(roots_u)
        return RootBunch(
            constraint = root_constraints.mirror_imag,
            c = roots_c,
            i = roots_i,
            u = roots_u,
        )

    def MI2MQ(self, roots_c, roots_i, roots_u):
        real_select = self.are_real(roots_c)
        roots_r = roots_c[real_select].real
        roots_c = roots_c[~real_select]

        pos_select = roots_c.imag > 0
        roots_c_neg = roots_c[~pos_select]
        roots_c_pos = roots_c[pos_select]
        rPB = TFmath.nearest_pairs(roots_c_pos, roots_c_neg.conjugate())
        roots_u2 = list(rPB.l1_remain) + [r.conjugate() for r in rPB.l2_remain]
        roots_c = []
        for r1, r2 in rPB.r12_list:
            if not self.strict or self.are_same(r1, r2):
                #roots_c.append((r1 + r2) / 2)
                #TODO, this seems to work better, not clear why..
                roots_c.append(r1)
            else:
                roots_u2.append(r1)
                roots_u2.append(r2.conjugate())
        roots_c = np.array(roots_c)
        roots_u = np.concatenate([roots_u, roots_u2])

        select_zero = self.are_zero(roots_i)
        roots_i = roots_i[~select_zero]
        z = np.count_nonzero(select_zero)
        return RootBunch(
            constraint = root_constraints.mirror_quad,
            z = z,
            c = roots_c,
            r = roots_r,
            i = roots_i,
            u = roots_u,
        )

    def MR2MQ(self, roots_c, roots_r, roots_u):
        imag_select = self.are_imag(roots_c)
        roots_i = roots_c[imag_select].imag
        roots_c = roots_c[~imag_select]

        pos_select = roots_c.real > 0
        roots_c_neg = roots_c[~pos_select]
        roots_c_pos = roots_c[pos_select]
        rPB = TFmath.nearest_pairs(roots_c_pos, -roots_c_neg.conjugate())
        roots_u2 = list(rPB.l1_remain) + [-r.conjugate() for r in rPB.l2_remain]
        roots_c = []
        for r1, r2 in rPB.r12_list:
            if not self.strict or self.are_same(r1, r2):
                #roots_c.append((r1 + r2) / 2)
                #TODO, this seems to work better, not clear why..
                roots_c.append(r1)
            else:
                roots_u2.append(r1)
                roots_u2.append(-r2.conjugate())
        roots_c = np.array(roots_c)
        roots_u = np.concatenate([roots_u, roots_u2])

        select_zero = self.are_zero(roots_r)
        roots_r = roots_r[~select_zero]
        z = np.count_nonzero(select_zero)
        return RootBunch(
            constraint = root_constraints.mirror_quad,
            z = z,
            c = roots_c,
            r = roots_r,
            i = roots_i,
            u = roots_u,
        )

