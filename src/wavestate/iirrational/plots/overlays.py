# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals

import numpy as np
import matplotlib.collections as mplcollect
#import matplotlib as mpl
#from matplotlib import patches
import matplotlib.transforms as mpltrans

from . import utilities

def plot_ZP_grab(
    self,
    fitter,
    duals,
    color = 'black',
    axB   = None,
):
    def rlog_F(r):
        F_Hz = r.imag
        BW = abs(r.real)
        #select = abs(r) < 1
        return BW, F_Hz

    ax = axB.mag2
    x_z = []
    y_z = []
    x_p = []
    y_p = []
    for b in duals:
        #coding_z = fitter.num_codings[b.idx_z]
        #z_rl, z_F = rlog_F(coding_z.roots_c(fitter)[0])
        z_rl, z_F = rlog_F(b.z)
        #coding_p = fitter.den_codings[b.idx_p]
        #p_rl, p_F = rlog_F(coding_p.roots_c(fitter)[0])
        p_rl, p_F = rlog_F(b.p)

        trans = mpltrans.composite_transform_factory(
            ax.transData,
            ax.transAxes.inverted(),
        )
        locs_z = trans.transform((z_F, z_rl))
        locs_p = trans.transform((p_F, p_rl))
        x_z.append(locs_z[0])
        y_z.append(locs_z[1])
        x_p.append(locs_p[0])
        y_p.append(locs_p[1])
    x_z = np.asarray(x_z)
    y_z = np.asarray(y_z)
    x_p = np.asarray(x_p)
    y_p = np.asarray(y_p)

    aspect = utilities.get_aspect(ax)
    r = .01 + ((x_z - x_p)**2 + (y_z - y_p)**2 / aspect**1)**.5,
    angles = np.angle((x_z - x_p) + 1j * (y_z - y_p) / aspect, deg = True)
    col = mplcollect.EllipseCollection(
        widths      = r,
        heights     = .012,
        angles      = angles,
        units       = 'width',
        facecolors  = 'none',
        edgecolors  = color,
        offsets     = np.vstack([(x_z + x_p)/2, (y_z + y_p)/2]).T,
        transOffset = ax.transAxes,
    )
    ax.add_collection(col)
    return axB

