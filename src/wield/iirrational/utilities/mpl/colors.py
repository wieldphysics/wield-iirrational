#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""


import os
from os import path

from wield import declarative
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# from http://godsnotwheregodsnot.blogspot.com/2012/09/color-distribution-methodology.html

color_array = (
    np.array(
        [
            [1, 0, 103],  # 0
            [158, 0, 142],  # 1
            [14, 76, 161],  # 2
            [0, 95, 57],  # 3
            [149, 0, 58],  # 4
            [255, 147, 126],  # 5
            [164, 36, 0],  # 6
            [0, 21, 68],  # 7
            [98, 14, 0],  # 8
            [0, 0, 255],  # 3
            [255, 0, 86],  # 2
            [0, 125, 181],  # 1
            [106, 130, 108],  # 1
            [0, 174, 126],  # 1
            [194, 140, 159],  # 1
            [190, 153, 112],  # 1
            [0, 143, 156],  # 1
            [95, 173, 78],  # 1
            [255, 0, 0],  # 1
            [255, 0, 246],  # 1
            [255, 2, 157],  # 1
            [104, 61, 59],  # 1
            [255, 116, 163],  # 1
            [150, 138, 232],  # 1
            [107, 104, 130],  # 9
            [152, 255, 82],  # 1
            [167, 87, 64],  # 1
            [1, 255, 254],  # 1
            [255, 238, 232],  # 1
            [254, 137, 0],  # 1
            [145, 208, 203],  # 9
            [189, 198, 255],  # 1
            [1, 208, 255],  # 1
            [187, 136, 0],  # 1
            [117, 68, 177],  # 1
            [165, 255, 210],  # 1
            [255, 166, 254],  # 1
            [119, 77, 0],  # 1
            [122, 71, 130],  # 1
            [0, 255, 0],  # 0
            [38, 52, 0],  # 1
            [0, 71, 84],  # 1
            [67, 0, 44],  # 1
            [181, 0, 255],  # 1
            [255, 177, 103],  # 1
            [255, 219, 102],  # 1
            [144, 251, 146],  # 1
            [213, 255, 0],  # 2
            [126, 45, 210],  # 1
            [189, 211, 147],  # 1
            [229, 111, 254],  # 1
            [222, 255, 116],  # 1
            [0, 255, 120],  # 1
            [0, 155, 255],  # 1
            [0, 100, 1],  # 1
            [255, 229, 2],  # 6
            [0, 118, 255],  # 1
            [133, 169, 0],  # 1
            [0, 185, 23],  # 1
            [120, 130, 49],  # 1
            [0, 255, 198],  # 1
            [255, 110, 65],  # 1
            [232, 94, 190],  # 1
        ]
    )
    / 255.0
)
