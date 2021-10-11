# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np


def SNR_fix(aid):
    if aid.hint('trust_SNR'):
        return
    with aid.log_heading('SNR Fix Test'):
        W = aid.fitter.W
        SNR_max = aid.hint('SNR_max')
        if SNR_max is not None:
            Wclip = np.minimum(W, SNR_max)
            if np.any(Wclip < W):
                W = Wclip
                aid.fitter.W = Wclip
                aid.invalidate_fitters()
                aid.log_warn(
                    3,
                    """
                    Applying 'SNR_max'={} to clip the weights
                    """.format(SNR_max),
                )

                #TODO, this validate=False is for relative_degree checks
                #since this is so early. Need a better mehtod
                aid.fitter_update(
                    representative = False,
                    validate = False,
                )

        effective_data_points = (np.sum(W**2))**2 / np.sum(W**4)
        ratio = effective_data_points / len(W)
        Wmax = max(W)

        pref_ratio = aid.hint('SNR_regularize_ratio')
        pref_scale = aid.hint('SNR_regularize_scale')

        did_warn = False
        if ratio < 1 - pref_scale / Wmax:
            aid.log_warn(
                3,
                """
                The number of effective data points N=(ΣW^2)^2/(ΣW^4)={:.2e}*len(W)
                [where W=SNR] is below the configured 'SNR_regularize_scale'={},
                given the maximum SNR={}. Now Finding an SNR ceiling that balances
                the ratio with max SNR.
                """.format(ratio, pref_scale, Wmax),
            )
            did_warn = True
        elif ratio < pref_ratio:
            aid.log_warn(
                3,
                """
                The number of effective data points N=(ΣW^2)^2/(ΣW^4)={:.2e}*len(W)
                [where W=SNR] is below the configured 'SNR_regularize_ratio'={}.
                Now Finding an SNR ceiling to adjust to that ratio.
                """.format(ratio, pref_ratio),
            )
            did_warn = True

        if did_warn:
            aid.log_rationale(
                3,
                """
                Disable this adjustment by using a smaller value of
                'SNR_regularize_ratio', setting 'fix_SNR'=False, or
                setting 'trust_SNR'=True.
                """.format(pref_ratio),
            )
        else:
            return

        W = np.sort(aid.fitter.W)
        idx_min = 0
        idx_max = len(W)

        #binary search
        while idx_min+1 < idx_max:
            idx_mid = (idx_max - idx_min) // 2 + idx_min
            W_ceiling = W[idx_mid]
            if W_ceiling <= 0:
                idx_min = idx_mid
                continue
            W_new = np.minimum(W_ceiling, W)
            effective_data_points = (np.sum(W_new**2))**2 / np.sum(W_new**4)
            ratio = effective_data_points / len(W)
            pref_scale_ratio = (1 - pref_scale / W_ceiling)
            if ratio < pref_ratio or ratio < pref_scale_ratio:
                idx_max = idx_mid
            else:
                idx_min = idx_mid
            aid.log_debug(
                9,
                'W:ratio',
                pref_ratio,
                pref_scale_ratio,
                ratio,
                idx_min,
                idx_max,
                W_ceiling
            )

        W_ceiling = W[idx_min]
        aid.log_warn(
            3,
            """
            Using SNR<{} ceiling.
            """.format(W_ceiling),
        )
        aid.fitter.W = np.minimum(aid.fitter.W, W_ceiling)
        aid.invalidate_fitters()

        #TODO, this validate=False is for relative_degree checks
        #since this is so early. Need a better mehtod
        aid.fitter_update(
            representative = False,
            validate = False,
        )


def check_sample_variance(aid):
    #TODO
    return

def check_sample_variance_adv(aid):
    #TODO
    return

