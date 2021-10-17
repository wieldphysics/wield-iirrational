# -*- coding: utf-8 -*-
"""
"""


from .. import representations
from .mappings import coding_maps

RBalgo = representations.RBalgo


def ZPKrep2MRF(
    ZPKrep,
    coding_map          = coding_maps.nlFBW_safe,
    delay_s             = None,
    insert_poles        = None,
    insert_zeros        = None,
    distance_limit_auto = None,
    **kwargs
):
    try:
        repdict = coding_map.rep
    except AttributeError:
        coding_submap = coding_map
    else:
        try:
            if ZPKrep.F_nyquist_Hz is None:
                rep = 'Sf'
            else:
                rep = 'Z'
            coding_submap = repdict[rep]
        except KeyError:
            raise RuntimeError((
                "Conversion to representation \"{}\" not supported"
                " in used coding_map."
            ).format(ZPKrep.rep))

    gain_coding = coding_submap.gain()
    gain_coding.setup(gain = ZPKrep.gain)

    if delay_s is None:
        delay_s = ZPKrep.delay_s

    #pregenerate the fitter to use as the system source for the codings
    fitterMRF = coding_submap.mrf_default(
        data                = ZPKrep.data,
        W                   = ZPKrep.W,
        F_Hz                = ZPKrep.F_Hz,
        zeros_overlay       = ZPKrep.zeros_overlay,
        poles_overlay       = ZPKrep.poles_overlay,
        F_nyquist_Hz        = ZPKrep.F_nyquist_Hz,
        residuals_log_im_scale = ZPKrep.residuals_log_im_scale,
        delay_s             = delay_s,
        num_codings         = [],
        den_codings         = [],
        gain_coding         = gain_coding,
        coding_map          = coding_submap,
        distance_limit_auto = distance_limit_auto,
        **kwargs
    )
    fitterMRF.delay_s = ZPKrep.delay_s

    #TODO, reset the gain to use the gain_effect parameters
    num_codings = []
    den_codings = []
    Z = RBalgo.expect(
        ZPKrep.zeros,
        constraint = representations.root_constraints.mirror_real
    )
    P = RBalgo.expect(
        ZPKrep.poles,
        constraint = representations.root_constraints.mirror_real
    )
    if insert_zeros is not None:
        #TODO, add completion?
        Z = Z * RBalgo.expect(
            insert_zeros,
            constraint = representations.root_constraints.mirror_real
        )
    if insert_poles is not None:
        #TODO, add completion?
        P = P * RBalgo.expect(
            insert_poles,
            constraint = representations.root_constraints.mirror_real
        )

    num_codings = []
    num_codings.extend(
        encode(
            fitterMRF,
            Z.r,
            coding_submap.num_r,
            coding_submap.num_r_u,
            coding_submap.num_collect_r,
            coding_submap.is_unstable,
        )
    )
    num_codings.extend(
        encode(
            fitterMRF,
            Z.c,
            coding_submap.num_c,
            coding_submap.num_c_u,
            coding_submap.num_collect_c,
            coding_submap.is_unstable,
        )
    )

    den_codings = []
    den_codings.extend(
        encode(
            fitterMRF,
            P.r,
            coding_submap.den_r,
            coding_submap.den_r_u,
            coding_submap.den_collect_r,
            coding_submap.is_unstable,
        )
    )
    den_codings.extend(
        encode(
            fitterMRF,
            P.c,
            coding_submap.den_c,
            coding_submap.den_c_u,
            coding_submap.den_collect_c,
            coding_submap.is_unstable,
        )
    )

    #now set the codings
    fitterMRF.num_codings = num_codings
    fitterMRF.den_codings = den_codings
    #set the gain again to account for gain_effect
    fitterMRF.gain = ZPKrep.gain
    #TODO, HACK
    return fitterMRF


def encode(fitter, r_list, st_t, us_t, collect, is_unstable):
    """
    Internal function to perform the encoding
    """
    if us_t is not None:
        r_list_st = []
        r_list_us = []
        for root in r_list:
            if is_unstable(root):
                r_list_us.append(root)
            else:
                r_list_st.append(root)
    else:
        r_list_st = list(r_list)
        r_list_us = []

    if us_t is None:
        us_t = st_t

    codings = []
    while r_list_st:
        coding = st_t(fitter)
        roots = r_list_st[-collect:]
        r_list_st[-collect:] = []
        coding.update_roots(*roots)
        codings.append(coding)
    while r_list_us:
        coding = us_t(fitter)
        roots = r_list_us[-collect:]
        r_list_us[-collect:] = []
        coding.update_roots(*roots)
        codings.append(coding)
    return codings


def MRF2MRF(
    fitter,
    ZPKrep     = None,
    zeros      = None,
    poles      = None,
    gain       = None,
    p_c        = None,
    z_c        = None,
    p_r        = None,
    z_r        = None,
    delay_s    = None,
    coding_map = None,
    **kwargs
):

    if ZPKrep is None:
        ZPKrep = fitter.ZPKrep

    if zeros is None:
        zeros = ZPKrep.zeros
    else:
        zeros = fitter.RBalgo.expect(zeros)

    if poles is None:
        poles = ZPKrep.poles
    else:
        poles = fitter.RBalgo.expect(poles)

    if z_r is None:
        z_r = zeros.r
    if p_r is None:
        p_r = poles.r
    if z_c is None:
        z_c = zeros.c
    if p_c is None:
        p_c = poles.c

    if gain is None:
        gain = ZPKrep.gain

    zeros = representations.asMRRB(r = z_r, c = z_c)
    poles = representations.asMRRB(r = p_r, c = p_c)
    coding_map, num_codings, den_codings = ZP2codings(
        fitter, zeros, poles,
        coding_map = coding_map
    )

    gain_coding = coding_map.gain()
    gain_coding.setup(gain = gain)

    if delay_s is None:
        delay_s = ZPKrep.delay_s

    #pregenerate the fitter to use as the system source for the codings
    fitterMRF = coding_map.mrf_default(
        parent        = fitter,
        data          = ZPKrep.data,
        F_Hz          = ZPKrep.F_Hz,
        W             = ZPKrep.W,
        num_codings   = num_codings,
        den_codings   = den_codings,
        gain_coding   = gain_coding,
        codings       = coding_map.module,
        zeros_overlay = ZPKrep.zeros_overlay,
        poles_overlay = ZPKrep.poles_overlay,
        delay_s       = delay_s,
        **kwargs
    )
    return fitterMRF


def ZP2codings(
    fitter, zeros, poles,
    coding_map = None,
    **kwargs
):
    if coding_map is None:
        coding_map = fitter.coding_map

    try:
        repdict = coding_map.rep
    except AttributeError:
        coding_map = coding_map
    else:
        if fitter.F_nyquist_Hz is None:
            rep = 'Sf'
        else:
            rep = 'Z'
        coding_map = repdict[rep]

    zeros = fitter.RBalgo.expect(zeros, fitter.root_constraint)
    poles = fitter.RBalgo.expect(poles, fitter.root_constraint)

    num_codings = []
    num_codings.extend(
        encode(
            fitter,
            zeros.r,
            coding_map.num_r,
            coding_map.num_r_u,
            coding_map.num_collect_r,
            coding_map.is_unstable,
        )
    )
    num_codings.extend(
        encode(
            fitter,
            zeros.c,
            coding_map.num_c,
            coding_map.num_c_u,
            coding_map.num_collect_c,
            coding_map.is_unstable,
        )
    )
    den_codings = []
    den_codings.extend(
        encode(
            fitter,
            poles.r,
            coding_map.den_r,
            coding_map.den_r_u,
            coding_map.den_collect_r,
            coding_map.is_unstable,
        )
    )
    den_codings.extend(
        encode(
            fitter,
            poles.c,
            coding_map.den_c,
            coding_map.den_c_u,
            coding_map.den_collect_c,
            coding_map.is_unstable,
        )
    )

    return coding_map, num_codings, den_codings
