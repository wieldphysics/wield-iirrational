"""
"""
from __future__ import division, print_function, unicode_literals

import declarative
from . import codings_s
from . import codings_z

_common = dict(
)

coding_maps = declarative.Bunch(
    RI = declarative.Bunch(
        rep = declarative.Bunch(
            Sf = codings_s.coding_maps.RI,
            Z  = codings_z.coding_maps.RI,
        ),
        **_common
    ),
    FBW = declarative.Bunch(
        rep = declarative.Bunch(
            Sf = codings_s.coding_maps.FBW,
            Z  = codings_z.coding_maps.FBW,
        ),
        **_common
    ),
    nlFBW = declarative.Bunch(
        rep = declarative.Bunch(
            Sf = codings_s.coding_maps.nlFBW,
            Z  = codings_z.coding_maps.nlFBW,
        ),
        **_common
    ),
    nlFBW_safe = declarative.Bunch(
        rep = declarative.Bunch(
            Sf = codings_s.coding_maps.nlFBW_safe,
            Z  = codings_z.coding_maps.nlFBW_safe,
        ),
        **_common
    ),
    SOS = declarative.Bunch(
        rep = declarative.Bunch(
            Sf = codings_s.coding_maps.SOS,
            Z  = codings_z.coding_maps.SOS,
        ),
        **_common
    ),
    #SOS_safe = declarative.Bunch(
    #    rep = declarative.Bunch(
    #        Sf = codings_s.coding_maps.SOS_safe,
    #        Z  = codings_z.coding_maps.SOS_safe,
    #    ),
    #    **_common
    #),
)
