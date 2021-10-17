"""
Contains setup functions for test data, these are returned in a QuickBunch dictionary, with some annotation about the number of data sets
"""

import numpy as np

import declarative
from ..utilities.np import logspaced
from . import utilities


def meta_generator(
        SNR = 10,
        log = False,
        N   = 1000,
        Nreal_max = 10,
        Ncplx_max = 10,
):
    def dataset_gen(
            sN  = 0,
            iN  = 0,
            **kwargs
    ):
        rand = np.random.RandomState()
        rand.seed(sN)

        Nreal = rand.randint(Nreal_max) + 1
        Ncplx = rand.randint(Ncplx_max) + 1

        NrealP = (Nreal + 1)//2  # rand.randint(Nreal//2 + 1) + (Nreal + 1)//4
        NrealZ = Nreal - NrealP
        NcplxP = (Ncplx + 1)//2  # (Ncplx + 1)//4 + rand.randint((Ncplx + 1)//2 + 1)
        NcplxZ = Ncplx - NcplxP

        F_max = 100
        F_max_R = 80
        F_min = .001
        F_nyquist_Hz = None
        if log:
            F_Hz = logspaced(F_min, F_max, N)
            def densityBW(F_Hz):
                return np.log(F_max / F_min) * F_Hz / N
            Preal = -F_min * np.exp(rand.uniform(0, 1, NrealP) * np.log(F_max_R / F_min))
            PcplxF = 2 * F_min * np.exp(rand.uniform(0, 1, NcplxP) * np.log(F_max_R / F_min / 2))
            PcplxR = -np.exp(rand.uniform(0, 1, NcplxP)) * densityBW(PcplxF)
            Pcplx = PcplxR + 1j * PcplxF

            Zreal = -F_min * np.exp(rand.uniform(0, 1, NrealZ) * np.log(F_max_R / F_min))
            ZcplxF = 2 * F_min * np.exp(rand.uniform(0, 1, NcplxZ) * np.log(F_max_R / F_min / 2))
            ZcplxR = -np.exp(rand.uniform(0, 2, NcplxZ)) * densityBW(ZcplxF)
            Zcplx = ZcplxR + 1j * ZcplxF
        else:
            F_Hz = np.linspace(0, F_max, N)
            def densityBW(F_Hz):
                return (F_max - F_min) / N
            Preal = -F_min * np.exp(rand.uniform(0, 1, NrealP) * np.log(F_max_R / F_min))
            Zreal = -F_min * np.exp(rand.uniform(0, 1, NrealZ) * np.log(F_max_R / F_min))
            PcplxF = rand.uniform(0, F_max_R, NcplxP)
            ZcplxF = rand.uniform(0, F_max_R, NcplxZ)
            PcplxR = -rand.uniform(0, densityBW(PcplxF), NcplxP)
            ZcplxR = -rand.uniform(0, densityBW(ZcplxF), NcplxZ)
            Pcplx = PcplxR + 1j * PcplxF
            Zcplx = ZcplxR + 1j * ZcplxF
        K = 1
        Z = np.concatenate([Zreal, Zcplx, Zcplx.conjugate()])
        P = np.concatenate([Preal, Pcplx, Pcplx.conjugate()])

        Rp = np.prod(abs(Pcplx)) * np.prod(abs(Preal))
        Rz = np.prod(abs(Zcplx)) / np.prod(abs(Zreal))
        R = Rp / Rz
        print("RATIO1", R, Rp, Rz)
        if abs(np.log10(Rp)) > 5:
            Pcplx /= Rp**(1./len(Pcplx))
        if abs(np.log10(Rz)) > 5:
            Zcplx /= Rz**(1./len(Zcplx))
        Rp = np.prod(abs(Pcplx)) * np.prod(abs(Preal))
        Rz = np.prod(abs(Zcplx)) / np.prod(abs(Zreal))
        R = Rp / Rz
        print("RATIO2", R, Rp, Rz)

        return utilities.generator_autofill(
            F_Hz          = F_Hz,
            SNR           = SNR,
            F_nyquist_Hz  = F_nyquist_Hz,
            bestfit_ZPK_s = [Z, P, K],
            sN            = sN,
            iN            = iN,
        )
    return dataset_gen


datasets = wavestate.bunch.Bunch()

datasets.rand10_log1k = utilities.make_description(
    generator = meta_generator(
        SNR = 10,
        log = True,
        N   = 1000,
        Nreal_max = 5,
        Ncplx_max = 10,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand10_lin1k = utilities.make_description(
    generator = meta_generator(
        SNR = 10,
        log = False,
        N   = 1000,
        Nreal_max = 5,
        Ncplx_max = 10,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand6_log100 = utilities.make_description(
    generator = meta_generator(
        SNR = 10,
        log = True,
        N   = 100,
        Nreal_max = 4,
        Ncplx_max = 6,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand6_lin100 = utilities.make_description(
    generator = meta_generator(
        SNR = 10,
        log = False,
        N   = 100,
        Nreal_max = 4,
        Ncplx_max = 6,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand2_log100 = utilities.make_description(
    generator = meta_generator(
        SNR = 10,
        log = True,
        N   = 100,
        Nreal_max = 2,
        Ncplx_max = 2,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand2_lin100 = utilities.make_description(
    generator = meta_generator(
        SNR = 10,
        log = False,
        N   = 100,
        Nreal_max = 2,
        Ncplx_max = 2,
    ),
    sets      = -1,
    instances = -1,
)

datasets.rand6_log100E = utilities.make_description(
    generator = meta_generator(
        SNR = None,
        log = True,
        N   = 100,
        Nreal_max = 4,
        Ncplx_max = 6,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand6_lin100E = utilities.make_description(
    generator = meta_generator(
        SNR = None,
        log = False,
        N   = 100,
        Nreal_max = 6,
        Ncplx_max = 6,
    ),
    sets      = -1,
    instances = -1,
)

datasets.rand2_log100E = utilities.make_description(
    generator = meta_generator(
        SNR = None,
        log = True,
        N   = 100,
        Nreal_max = 2,
        Ncplx_max = 2,
    ),
    sets      = -1,
    instances = -1,
)


datasets.rand2_lin100E = utilities.make_description(
    generator = meta_generator(
        SNR = None,
        log = False,
        N   = 100,
        Nreal_max = 2,
        Ncplx_max = 2,
    ),
    sets      = -1,
    instances = -1,
)
