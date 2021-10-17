"""
Contains setup functions for test data, these are returned in a QuickBunch dictionary, with some annotation about the number of data sets
"""

import numpy as np
import declarative

from ..utilities.np import logspaced
from . import utilities


def meta_generator(
        Z, P, K,
        N            = 1000,
        SNR_default  = 10,
        delay_s      = 0,
        F_nyquist_Hz = None,
):
    def dataset_gen(
            sN           = 0,
            iN           = 0,
            SNR          = None,
            invert       = False,
            F_nyquist_Hz = F_nyquist_Hz,
            N            = N,
            **kwargs
    ):
        if SNR is None:
            SNR = SNR_default

        if not invert:
            ZPK_s = (Z, P, K)
        else:
            ZPK_s = (P, Z, 1/K)

        if sN == 0:
            F_Hz = np.linspace(0, 10, N)
        elif sN == 1:
            F_Hz = logspaced(.01, 10, N)
        elif sN == 2:
            F_Hz = np.linspace(0, 100, N)
        else:
            F_Hz = logspaced(.01, 100, N)

        return utilities.generator_autofill(
            F_Hz          = F_Hz,
            SNR           = SNR,
            F_nyquist_Hz  = F_nyquist_Hz,
            bestfit_ZPK_s = ZPK_s,
            delay_s       = delay_s,
            sN            = sN,
            iN            = iN,
        )
    return dataset_gen


datasets = wavestate.bunch.Bunch()

datasets.simple0 = utilities.make_description(
    generator = meta_generator(
        (),
        (-1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple0E = utilities.make_description(
    generator = meta_generator(
        (),
        (-1,),
        1,
        N = 100,
        SNR_default = None,
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple1 = utilities.make_description(
    generator = meta_generator(
        (-1, -1, -1, -1, -1),
        (-.1, -.1, -.1),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple1E = utilities.make_description(
    generator = meta_generator(
        (-1, -1, -1, -1, -1),
        (-.1, -.1, -.1),
        1,
        N = 30,
        SNR_default = None,
    ),
    sets      = 4,
    instances = -1,
)

datasets.delay = utilities.make_description(
    generator = meta_generator(
        (), (), 1, delay_s = .01,
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple2 = utilities.make_description(
    generator = meta_generator(
        (-1, -1, -.1 + 1.1j, -.1 - 1.1j, -1),
        (-.1, -.1, -.1 + 3j, -.1 - 3j, -.01 + 1j, -.01 - 1j, -1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple2_us = utilities.make_description(
    generator = meta_generator(
        (-1, -1, -.1 + 1.1j, -.1 - 1.1j, -1),
        (-.1, -.1, +.1 + 3j, +.1 - 3j, -.01 + 1j, -.01 - 1j, -1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple2_nmd = utilities.make_description(
    generator = meta_generator(
        (-1, -1, +.1 + 1.1j, +.1 - 1.1j, -1),
        (-.1, -.1, -.1 + 3j, -.1 - 3j, -.01 + 1j, -.01 - 1j, -1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple2_us2 = utilities.make_description(
    generator = meta_generator(
        (-1, -1, -.1 + 1.1j, -.1 - 1.1j, -1),
        (-.1, -.1, +.1 + 3j, +.1 - 3j, +.01 + 1j, +.01 - 1j, -1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple2lownyquist = utilities.make_description(
    generator = meta_generator(
        (-1, -1, -.1 + 1.1j, -.1 - 1.1j, -1),
        (-.1, -.1, -.1 + 3j, -.1 - 3j, -.01 + 1j, -.01 - 1j, -1,),
        1,
        F_nyquist_Hz = 10,
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple3 = utilities.make_description(
    generator = meta_generator(
        (-2, -2,),
        (-.1,  -.1, -.1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.simple4 = utilities.make_description(
    generator = meta_generator(
        (-.1,),
        (-1, -1,),
        1
    ),
    sets      = 4,
    instances = -1,
)

datasets.testZPK4x = utilities.make_description(
    generator = meta_generator(
        (+.1, +1 + 1j, +1 - 1j),
        (-.5, -2 + 1j, -2 - 1j, -2 + 2j, -2 - 2j),
        1
    ),
    sets      = 1,
    instances = -1,
)

datasets.testZPK4xHQ = utilities.make_description(
    #these are tuned to be right over 2x the log spacing
    #which is a special spacing used in the MRF code to not modify root locations
    generator = meta_generator(
        (),  # (-.1, -.1),
        (-.03 + 2.0j, -.03 - 2.0j),
        1,
        SNR_default = None,
    ),
    sets      = 1,
    instances = -1,
)
