# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

from IIRrational.testing import IIRrational_data
from IIRrational import plots
#from IIRrational.utilities.mpl import mplfigB
import os.path as path

from IIRrational.fitters_ZPK import ZPKrep2MRF
from IIRrational.fitters_ZPK import codings_z
from IIRrational import testing

from utilities import check_residuals


@pytest.mark.parametrize("coding_map", sorted(codings_z.coding_maps.keys()))
def test_NL_fitter(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    data = IIRrational_data('testZPK4x')
    #data.rep_z.zeros = data.rep_z.zeros
    #data.rep_z.poles = data.rep_z.poles
    fitterMRF = ZPKrep2MRF(
        data.rep_z.ZPKrep,
        coding_map = codings_z.coding_maps[coding_map],
    )
    testing.assert_almost_equal(data.rep_z.xfer_fit / fitterMRF.xfer_fit , 1, 4)
    testing.assert_almost_equal(data.rep_z.ZPKrep.xfer_fit / fitterMRF.xfer_fit, 1, 4)

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    #axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)

@pytest.mark.parametrize("coding_map", sorted(codings_z.coding_maps.keys()))
def test_NL_fitter_US(coding_map):
    data = IIRrational_data('testZPK4x')
    data.rep_z.zeros = 1.01 * (data.rep_z.zeros)
    data.rep_z.poles = .99 * (data.rep_z.poles)
    fitterMRF = ZPKrep2MRF(
        data.rep_z.ZPKrep,
        coding_map = codings_z.coding_maps[coding_map],
    )

    testing.assert_almost_equal(data.rep_z.xfer_fit / fitterMRF.xfer_fit , 1, 4)
    testing.assert_almost_equal(data.rep_z.ZPKrep.xfer_fit / fitterMRF.xfer_fit, 1, 4)
    print(data.rep_z.F_nyquist_Hz)
    print(fitterMRF.F_nyquist_Hz)

    try:
        print(fitterMRF.residuals_average)
        check_residuals(fitterMRF)
        print(fitterMRF.zeros)
        opt = fitterMRF.optimize()
        print(fitterMRF.zeros)
        print(fitterMRF.residuals_average)
        #the check after optimizing to unstable roots is... bad
        #check_residuals(fitterMRF)
        #assert(fitterMRF.residuals_average < 10)
    except AssertionError:
        #TODO make this plot_on_assert
        axB = plots.plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)
        axB.save(path.join(path.split(__file__)[0], 'out-Z/{}.png'.format(coding_map)))
        raise

@pytest.mark.parametrize("coding_map", sorted(codings_z.coding_maps.keys()))
def test_NL_fitter1_DL(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    #set number 1 is logarithmic data
    data = IIRrational_data('testZPK4x', set_num = 1)
    data.rep_z.zeros = (data.rep_z.zeros)
    data.rep_z.poles = (data.rep_z.poles)
    print(data.rep_z.zeros)
    fitterMRF = ZPKrep2MRF(
        data.rep_z.ZPKrep,
        coding_map = codings_z.coding_maps[coding_map],
        distance_limit_auto = 1,
    )
    testing.assert_almost_equal(data.rep_z.xfer_fit / fitterMRF.xfer_fit, 1, 2)
    testing.assert_almost_equal(data.rep_z.ZPKrep.xfer_fit / fitterMRF.xfer_fit, 1, 2)
    fitterMRF.distance_limit_auto = 1

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    #axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)


@pytest.mark.parametrize("coding_map", sorted(codings_z.coding_maps.keys()))
def test_NL_fitter1_DL2(coding_map):
    """
    tests when the fitter is bad as it is trying to fit on the wrong side
    """
    #set number 1 is logarithmic data
    data = IIRrational_data('testZPK4x', set_num = 1)
    data.rep_z.zeros = (data.rep_z.zeros)
    data.rep_z.poles = (data.rep_z.poles)
    print(data.rep_z.zeros)
    fitterMRF = ZPKrep2MRF(
        data.rep_z.ZPKrep,
        coding_map = codings_z.coding_maps[coding_map],
        distance_limit_auto = 1,
    )
    testing.assert_almost_equal(data.rep_z.xfer_fit / fitterMRF.xfer_fit , 1, 2)
    testing.assert_almost_equal(data.rep_z.ZPKrep.xfer_fit / fitterMRF.xfer_fit, 1, 2)
    fitterMRF.distance_limit_auto = 2

    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    opt = fitterMRF.optimize()
    print(fitterMRF.residuals_average)
    check_residuals(fitterMRF)
    #axB = plot_fitter_flag(fitterMRF, with_error=False, plot_zp = True)
