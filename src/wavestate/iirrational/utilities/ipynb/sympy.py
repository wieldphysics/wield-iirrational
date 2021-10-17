# -*- coding: utf-8 -*-
"""
"""

import sympy
import phasor.math.dispatch_sympy
print("Sympy version: ", sympy.__version__)

#this makes the notebook sexy
sympy.init_printing(use_latex='mathjax')


from IPython.display import (
    display,
    display_pretty,
    display_html,
    display_jpeg,
    display_png,
    display_json,
    display_latex,
    display_svg,
    clear_output,
    Latex,
)


