#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


from . import dense


class StateSpaceBuilder(object):
    t_dense = dense.StateSpaceDense

    def ZPKdict(
        self,
        name,
        zdict,
        pdict,
        k,
        convention = 'scipy',
    ):
        ABCDEs = dense.zpkdict_cascade(
            zdict = zdict,
            pdict = pdict,
            k = k,
            convention = convention
        )
        SSs = []
        for idx, (A, B, C, D, E) in enumerate(ABCDEs):
            ss = self.t_dense(
                A = A,
                E = E,
                B = B,
                C = C,
                D = D,
                name = "TF_{}".format(idx)
            )
            SSs.append(ss)
        new = self.t_dense.chain(name, SSs)

        new.names_change('inputs', fr = "TF_0.i0".format(0), to = "{}.i0".format(name))
        new.names_change('output', fr = "TF_{}.o0".format(idx), to = "{}.o0".format(name))
        new.names_collect('states', to = name)
        new.names_collect('constr', to = name)
        return new

    def delay(self, name, delay_s, order = 1, method = 'bessel', **kwargs):
        if method == 'pade':
            A, B, C, D, E = dense.pade_delay(
                delay_s = delay_s,
                order = order,
                **kwargs
            )
            return self.t_dense(
                A = A, E = E, B = B, C = C, D = D,
                name = name,
            )
        if method == 'bessel':
            zdict, pdict, k = dense.bessel_delay(delay_s, order = order, **kwargs)
            return self.ZPKdict(name, zdict, pdict, k, **kwargs)
        else:
            raise RuntimeError("Unrecognized delay synthesis method")

        return
