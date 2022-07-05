# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from wavestate.iirrational import v2
from wavestate.iirrational.v2 import testing
from wavestate.iirrational import plots


try:
    from IIRrational_test_data import IIRrational_data_dev
except ImportError:
    module_import_skip = True
else:
    module_import_skip = False


@pytest.mark.skipif(module_import_skip, reason="cannot import IIRrational_test_data")
def test_quad_M0_P2V_10Hz_S(request, browser, plotsections, plot_verbosity):
    data_name = "quad_M0_P2V_10Hz"

    dataB = IIRrational_data_dev(
        data_name,
        instance_num=1,
        set_num=2,
    )

    out = v2.data2filter(
        F_Hz=dataB.F_Hz,
        data=dataB.data,
        SNR=dataB.SNR,
        F_nyquist_Hz=None,
        # order_initial = 80,
        total_degree_min=None,
        delay_s=0,
        mode="rational",
        baseline_only=True,
        # SNR_cuttoff = 50,
        hints=[
            testing.validate_plot_log(__file__, request),
        ],
    )

    # out.fitter.matched_pairs_clear(Q_rank_cutoff = .4)
    def plot():
        return plots.plot_fit(
            out.fitter,
            plot_zp=True,
            # plot_past_data = False,
            plot_past_data=True,
            xscale="log",
        )

    with testing.plot_on_assert(__file__, request, plot, plot_anyway=True):
        pass
        # assert(False)
    # with plot_on_assert(__file__, request, out.fitter, plot_anyway = True):
    #    pass
    return


# @pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
# def test_quad_M0_P2V_10Hz_rat_Z(browser, plotsections, plot_verbosity):
#    data_name = 'quad_M0_P2V_10Hz'
#    out = v1.rational_disc_fit(
#        IIRrational_data_dev(
#            data_name,
#            instance_num = 1,
#            set_num = 2,
#        ),
#        nyquist_final_Hz = 1e3,
#        order = 100,
#    )
#
#    import matplotlib.pyplot as plt
#    from wavestate.iirrational import plots
#    axB = plots.plot_fitter_flag(
#        out,
#        #don't use jpg as matplotlib can't write it out of the box (needs pillow)
#        fname = path.join(path.split(__file__)[0], 'out-cplx/test-rational-Z.png'),
#    )
#    plt.close(axB.fig)
#    plt.show()
#    return
#
# @pytest.mark.skipif(module_import_skip, reason = 'cannot import IIRrational_test_data')
# def test_quad_M0_P2V_10Hz_rat_S(browser, plotsections, plot_verbosity):
#    data_name = 'quad_M0_P2V_10Hz'
#    out = v1.rational_disc_fit(
#        IIRrational_data_dev(
#            data_name,
#            instance_num = 1,
#            set_num = 2,
#        ),
#        nyquist_final_Hz = 0,
#        order = 100,
#    )
#
#    import matplotlib.pyplot as plt
#    from wavestate.iirrational import plots
#    axB = plots.plot_fitter_flag(
#        out,
#        #don't use jpg as matplotlib can't write it out of the box (needs pillow)
#        fname = path.join(path.split(__file__)[0], 'out-cplx/test-rational-S.png'),
#    )
#    plt.close(axB.fig)
#    plt.show()
#    return
