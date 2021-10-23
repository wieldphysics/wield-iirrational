# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

from IIRrational.testing import IIRrational_data
from IIRrational import v1

try:
    import matlab.engine
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False

# @pytest.mark.skipif(module_import_skip, reason = 'cannot import matlab.engine')
@pytest.mark.skip(reason="Matlab testing not quite automated")
def test_simple1(browser, plotsections, plot_verbosity):
    mlab = matlab.engine.start_matlab("-desktop")
    data_name = "simple2"
    out = v1.data2filter(
        IIRrational_data(
            data_name,
            instance_num=1,
            set_num=2,
        ),
        order_initial=30,
        delay_s=0,
    )

    import matplotlib.pyplot as plt
    from IIRrational import plots

    axB = plots.plot_fitter_flag(
        out.fitter,
        # don't use jpg as matplotlib can't write it out of the box (needs pillow)
        fname=path.join(path.split(__file__)[0], "out-simple/test.png"),
    )
    plt.close(axB.fig)
    plt.show()
    return


if __name__ == "__main__":
    test_simple1(False, None, 10)
