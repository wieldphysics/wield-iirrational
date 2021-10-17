#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
Utilities to manipulate ZPK roots S to/from Z, and make transfer functions
"""


import numpy as np

from ..TFmath import match_SOS_pairs
from .. import representations


def filter2sortedZPK(filt):
    zpk = representations.asZPKTF(filt)
    z, p, k = zpk

    zzpp_pairs = match_SOS_pairs(
        zpk.zeros.r, zpk.zeros.c,
        zpk.poles.r, zpk.poles.c,
        F_nyquist_Hz = None
    )

    poles = []
    zeros = []

    for z1, z2, p1, p2 in zzpp_pairs:
        if z1 is not None:
            zeros.append(z1)
        if z2 is not None:
            zeros.append(z2)
        if p1 is not None:
            poles.append(p1)
        if p2 is not None:
            poles.append(p2)
    return zeros, poles, k

def filter2fotonZPK(
        filt,
        annotate_pairs = True,
        scale_gain     = True,
        plane          = None,
        zpk_output     = False,
):
    if plane is None:
        plane = 'f'
    zpk = representations.asZPKTF(filt)

    zzpp_pairs = match_SOS_pairs(
        zpk.zeros.r, zpk.zeros.c,
        zpk.poles.r, zpk.poles.c,
        F_nyquist_Hz = None
    )
    k = zpk.gain

    if scale_gain:
        if plane == 'f':
            k = k / ((np.pi*2)**(len(zpk.zeros) - len(zpk.poles)))
        elif plane == 'w':
            pass
        elif plane == 'n':
            def root_gain(r1, r2):
                k_adj = 1
                if r1 is not None:
                    r1 = -r1.conjugate()
                    aval = abs(r1)
                    if aval > 1e-9:
                        k_adj = aval

                if r2 is not None:
                    aval = abs(r2)
                    if aval > 1e-9:
                        k_adj = k_adj * aval
                    r2 = -r2.conjugate()

                return r1, r2, k_adj

            zzpp_pairs_mod = []
            for idx, (z1, z2, p1, p2) in enumerate(zzpp_pairs):
                z1, z2, k_adj = root_gain(z1, z2)
                k *= k_adj
                p1, p2, k_adj = root_gain(p1, p2)
                k /= k_adj
                zzpp_pairs_mod.append((z1, z2, p1, p2))
            zzpp_pairs = zzpp_pairs_mod

    if zpk_output:
        zero_list = []
        pole_list = []
        def root_collect(r1, r2):
            ret = []
            if r1 is not None:
                ret = [r1]
            else:
                ret = []

            if r2 is not None:
                ret += [r2]
            return ret

        for idx, (z1, z2, p1, p2) in enumerate(zzpp_pairs):
            zero_list.extend(root_collect(z1, z2,))
            pole_list.extend(root_collect(p1, p2,))
        return zero_list, pole_list, k

    def matlab_complex_str(num):
        if num is None:
            return None
        if(np.imag(num) == 0.):
            return ("{0}".format(num))
        else:
            r, i = np.real(num), np.imag(num)
            if i > 0:
                return ("{0} + {1}*i".format(r, i))
            else:
                return ("{0} - {1}*i".format(r, -i))

    if annotate_pairs:
        def root_strs(r1, r2, idx):
            s1 = matlab_complex_str(r1)
            s2 = matlab_complex_str(r2)
            if r2 is not None and r1 is not None:
                if len(s2) + len(s1) < 60:
                    return [
                        "{}; {};  % SOS {}AB".format(s1, s2, idx),
                    ]
                else:
                    return [
                        "{}; % SOS {}A".format(s1, idx),
                        "{}; % SOS {}B".format(s2, idx),
                    ]
            elif r1 is not None:
                return [
                    "{}; % SOS {}".format(s1, idx),
                ]
            elif r2 is not None:
                return [
                    "{}; % SOS {}".format(s2, idx),
                ]
            else:
                return []

    else:
        def root_strs(r1, r2, idx):
            s1 = matlab_complex_str(r1)
            s2 = matlab_complex_str(r2)
            if r2 is not None and r1 is not None:
                if len(s2) + len(s1) < 180:
                    return [
                        "{}; {};".format(s1, s2, idx),
                    ]
                else:
                    return [
                        "{};".format(s1, idx),
                        "{};".format(s2, idx),
                    ]
            elif r1 is not None:
                return [
                    "{};".format(s1, idx),
                ]
            elif r2 is not None:
                return [
                    "{};".format(s2, idx),
                ]
            else:
                return []

    pole_list = []
    zero_list = []
    for idx, (z1, z2, p1, p2) in enumerate(zzpp_pairs):
        zero_list.extend(root_strs(z1, z2, idx + 1,))
        pole_list.extend(root_strs(p1, p2, idx + 1,))

    ZPK_template = 'ZPK([{zeros}],[{poles}], {gain}, "{plane}")'
    if len(zero_list) == 0:
        zstr = ''
    else:
        zstr = '\n  ' + '\n  '.join(zero_list) + '\n'

    if len(pole_list) == 0:
        pstr = ''
    else:
        pstr = '\n  ' + '\n  '.join(pole_list) + '\n'

    return ZPK_template.format(
        zeros = zstr,
        poles = pstr,
        gain  = k,
        plane = plane,
    )


'''
def filter2matlabZPK(cascade):
    def matlab_complex_str(num):
        if(np.imag(num) == 0.):
            return ("{0}".format(num))
        else:
            return ("{0} + {1}*i".format(np.real(num), np.imag(num)))

    cascade = cascade.transformed_nyquist(F_nyquist=None)
    ZPK_template = 'zpk(2*pi*[;...\n{zeros}],2*pi*[;...\n{poles}],...\n{gain})'
    pole_list = []
    zero_list = []
    g_fix = 1
    zpk_num_mismatch = 0

    for pole in cascade.poles:
        if pole is not None:
            pole_list.append(matlab_complex_str(pole))
            g_fix *= pole
            zpk_num_mismatch += 1

    for zero in cascade.zeros:
        if zero is not None:
            zero_list.append(matlab_complex_str(zero))
            g_fix *= zero
            zpk_num_mismatch += 1

    return ZPK_template.format(
        zeros=';...\n\t' + ';...\n\t'.join(zero_list) + ';...\n',
        poles=';...\n\t' + ';...\n\t'.join(pole_list) + ';...\n',
        gain=cascade.gain * abs(g_fix) * (2*np.pi)**zpk_num_mismatch
    )


def filter2fotonSOS(cascade):
    raise NotImplementedError()
    if cascade.F_nyquist is None:
        print("Warning: cascade needs to be converted to Z domain for foton sos representation")
    SOS_template = 'sos({gain}, [{coeffs}])'
    coeff_list = []
    for idx in range(cascade.num_sos):
        coeff_list.append(cascade.coeff('b1', idx, default = True))
        coeff_list.append(cascade.coeff('b2', idx, default = True))
        coeff_list.append(cascade.coeff('a1', idx, default = True))
        coeff_list.append(cascade.coeff('a2', idx, default = True))
    def none_to_0(f):
        if f is None:
            return 0
        return f
    coeff_list = [str(none_to_0(i)) for i in coeff_list]
    return SOS_template.format(gain=cascade.gain, coeffs=';'.join(coeff_list))


def filter2fotonZRoots(filter):
    """

    """
    raise NotImplementedError()
    def matlab_complex_str(num):
        if(np.imag(num) == 0.):
            return ("{0}".format(num))
        else:
            return ("{0} - {1}*i".format(np.real(num), np.imag(num)))

    cascade = cascade.transformed_nyquist(F_nyquist=None)
    ZPK_template = 'ZPK([{zeros}],[{poles}],\n{gain}, "f")'
    pole_list = []
    zero_list = []
    for pole in cascade.poles:
        pole = -pole.conjugate()
        pole_list.append(matlab_complex_str(pole))
    for zero in cascade.poles:
        zero = -zero.conjugate()
        pole_list.append(matlab_complex_str(zero))

    return ZPK_template.format(
        zeros='\n\t' + ';\n\t'.join(zero_list) + '\n',
        poles='\n\t' + ';\n\t'.join(pole_list) + '\n',
        gain=cascade.gain
    )
'''
