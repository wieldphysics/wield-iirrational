# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

from IIRrational.testing import IIRrational_data
from IIRrational import v1


def test_simple1(browser, plotsections, plot_verbosity):
    #this opens a live plot, which we do not want, should incorporate that
    #functionality into plot_on_assert
    return
    data_name = 'simple2'
    out = v1.data2filter(
        IIRrational_data(
            data_name,
            instance_num = 1,
            set_num = 2,
        ),
        order_initial = 10,
        delay_s = 0,
    )

    import matplotlib.pyplot as plt
    from IIRrational import plots
    axB = plots.plot_fitter_flag(
        out.fitter,
        #don't use jpg as matplotlib can't write it out of the box (needs pillow)
        fname = path.join(path.split(__file__)[0], 'out-simple/test.png'),
    )
    plt.close(axB.fig)
    plt.show()
    return


if __name__=='__main__':
    test_simple1(False, None, 10)
