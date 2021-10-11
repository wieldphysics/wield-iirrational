"""
Some utilities to help with dataset creation and annotation
"""
from __future__ import division, print_function, unicode_literals
import declarative


import numpy as np

from .. import annotate
from .. import TFmath
from .. import representations
from .. import fitters_ZPK

from . import plots


def make_description(
    generator,
    sets        = 1,
    instances   = 1,
    description = None,
    annotation  = None,
    **kwargs
):
    doc = generator.__doc__
    if description is None:
        if doc is not None:
            description = annotate.padding_remove(doc)
        else:
            description = "<no description>"
    else:
        description = annotate.padding_remove(description)
        if annotation is None and doc is not None:
            annotation  = annotate.padding_remove(doc)
    return declarative.Bunch(
        generator   = generator,
        instances   = instances,
        description = description,
        annotation  = annotation,
        **kwargs
    )


def generator_autofill(
    F_Hz,
    SNR,
    F_nyquist_Hz,
    data               = None,
    bestfit_ZPK_z      = None,
    bestfit_ZPK_s      = None,
    delay_s            = 0,
    residuals_log_best = None,
    sN                 = 0,
    iN                 = 0,
    **kwargs
):
    if bestfit_ZPK_s:
        #replace with the ZPKs filled out fully
        bestfit_ZPK_s = TFmath.ZPK_fill(bestfit_ZPK_s)

    if bestfit_ZPK_z:
        #replace with the ZPKs filled out fully
        bestfit_ZPK_z = TFmath.ZPK_fill(bestfit_ZPK_z)

    if not bestfit_ZPK_z and F_nyquist_Hz is not None:
        if bestfit_ZPK_s:
            bestfit_ZPK_z = TFmath.StoZ(
                bestfit_ZPK_s,
                F_nyquist_Hz = F_nyquist_Hz,
            )

    if not bestfit_ZPK_s and F_nyquist_Hz is not None:
        if bestfit_ZPK_z:
            bestfit_ZPK_s = TFmath.ZtoS(
                bestfit_ZPK_z,
                F_nyquist_Hz = F_nyquist_Hz,
            )

    if data is None:
        if SNR is not None:
            rand = np.random.RandomState()
            rand.seed(iN)
            N = len(F_Hz)
            rel_noise = (rand.normal(1, 1/SNR, N) + 1j * rand.normal(0, 1/SNR, N))
        else:
            rel_noise = 1

        if bestfit_ZPK_s:
            data = rel_noise * TFmath.TF_ZPK(
                F_Hz,
                ZPK = bestfit_ZPK_s,
            )
        elif bestfit_ZPK_z:
            data = rel_noise * TFmath.TF_ZPK(
                F_Hz,
                ZPK = bestfit_ZPK_s,
                F_nyquist_Hz = F_nyquist_Hz,
            )

        if delay_s is not None:
            data = data * np.exp(-2j * np.pi * delay_s * F_Hz)

    if bestfit_ZPK_z is not None:
        rep_z = representations.ZPKwData(
            data         = data,
            F_Hz         = F_Hz,
            W            = SNR,
            F_nyquist_Hz = F_nyquist_Hz,
            ZPK          = bestfit_ZPK_z,
            delay_s      = delay_s,
        )
    else:
        rep_z = None

    if bestfit_ZPK_s is not None:
        rep_s = representations.ZPKwData(
            data         = data,
            F_Hz         = F_Hz,
            W            = SNR,
            F_nyquist_Hz = None,
            ZPK          = bestfit_ZPK_s,
            delay_s      = delay_s,
        )
    else:
        rep_s = None

    return declarative.Bunch(
        F_Hz               = F_Hz,
        data               = data,
        SNR                = SNR,
        F_nyquist_Hz       = F_nyquist_Hz,
        residuals_log_best = residuals_log_best,
        bestfit_ZPK_z      = bestfit_ZPK_z,
        bestfit_ZPK_s      = bestfit_ZPK_s,
        rep_z              = rep_z,
        rep_s              = rep_s,
        iN                 = iN,
        sN                 = sN,
        **kwargs
    )


def assert_almost_equal(arr1, arr2, decimals):
    np.testing.assert_allclose(arr1, arr2, decimals)
    #np.testing.assert_allclose(arr1.real, arr2.real, decimals)
    #np.testing.assert_allclose(arr1.imag, arr2.imag, decimals)


def sign_validate(aid, fitter):
    """
    To be added as a hint to data2filter
    """
    rep = fitter.ZPKrep
    xfer = rep.xfer_fit
    data = rep.data
    rat = data / xfer
    rat_ang = np.exp(1j * np.angle(rat))
    ang_avg_rep = np.sum(rat_ang * rep.W**2) / np.sum(rep.W**2)

    xfer = fitter.xfer_fit
    data = fitter.data
    rat = data / xfer
    rat_ang = np.exp(1j * np.angle(rat))
    ang_avg_fit = np.sum(rat_ang * fitter.W**2) / np.sum(fitter.W**2)
    #print("SGN: ", ang_avg_rep, ang_avg_fit)
    #axB = plots.plots.plot_fit(
    #    fitter,
    #    fname = 'test1.png',
    #)
    #axB = plots.plots.plot_fitter_flag(
    #    fitter,
    #    fname = 'test3.png',
    #)
    #axB = plots.plots.plot_fit(
    #    fitter.ZPKrep,
    #    fname = 'test2.png',
    #)

    if isinstance(fitter, fitters_ZPK.MultiReprFilterBase):
        for coding in list(fitter.num_codings) + list(fitter.den_codings):
            rB = representations.RootBunch(u = coding.roots(), constraint = representations.root_constraints.no_constraint)
            h1 = coding.transfer()
            h, lnG = rB.val_lnG(fitter.Xex_grid)
            h = h * np.exp(lnG)
            assert_almost_equal(h / h1, 1, 4)

    assert(ang_avg_fit.real > 0 and ang_avg_rep.real > 0)

sign_validate_hint = {
    'fitter_update_validate' : sign_validate,
    'fitter_check_validate' : sign_validate,
}

def rational_fitter_validate(rat_fitter, fitter):
    """
    To be added as a hint to data2filter
    """
    rep1 = rat_fitter.ZPKrep.xfer_fit
    rep2 = fitter.ZPKrep.xfer_fit
    rat = rat_fitter.xfer_fit
    mrf = fitter.xfer_fit

    assert_almost_equal(rep1 / rep2, 1, 5)
    assert_almost_equal(rat / rep2, 1, 5)
    assert_almost_equal(mrf / rep2, 1, 5)
    print("Checking Rational Fitter")

def sign_validate_and_plot_hint(pyfile, request):
    def sign_validate_plot(aid, fitter):
        with plots.plot_on_assert(pyfile, request, fitter, plot_anyway = False):
            try:
                sign_validate(aid, fitter)
            except AssertionError:
                print(fitter)
                print(fitter.F_nyquist_Hz)
                print(fitter.zeros)
                print(fitter.zeros_overlay)
                print(fitter.ZPKrep)
                #assert(False)
                raise
    hint = {
        'fitter_update_validate'   : sign_validate_plot,
        'fitter_check_validate'    : sign_validate_plot,
        'rational_fitter_validate' : rational_fitter_validate,
    }
    return hint

def stability_validate_and_plot_hint(pyfile, request):
    def sign_validate_plot(aid, fitter):
        with plots.plot_on_assert(pyfile, request, fitter, plot_anyway = False):
            assert(np.all(fitter.ZPKrep.poles.c.real <= 0))
            assert(np.all(fitter.ZPKrep.poles.r.real <= 0))
            try:
                sign_validate(aid, fitter)
            except AssertionError:
                print(fitter)
                print(fitter.F_nyquist_Hz)
                print(fitter.zeros)
                print(fitter.zeros_overlay)
                print(fitter.ZPKrep)
                #assert(False)
                raise
    hint = {
        'fitter_update_validate'   : sign_validate_plot,
        'fitter_check_validate'    : sign_validate_plot,
        'rational_fitter_validate' : rational_fitter_validate,
    }
    return hint

def sign_validate_and_digest_hint(pyfile, request):
    def sign_validate_plot(aid, fitter):
        with plots.digest_on_assert(pyfile, request, aid, plot_anyway = False):
            try:
                sign_validate(aid, fitter)
            except AssertionError:
                print(fitter)
                print(fitter.F_nyquist_Hz)
                print(fitter.zeros)
                print(fitter.zeros_overlay)
                print(fitter.ZPKrep)
                #assert(False)
                raise
    hint = {
        'fitter_update_validate'   : sign_validate_plot,
        'fitter_check_validate'    : sign_validate_plot,
        'rational_fitter_validate' : rational_fitter_validate,
    }
    return hint





