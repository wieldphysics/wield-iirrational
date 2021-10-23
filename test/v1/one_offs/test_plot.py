# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import numpy as np
import os.path as path
import matplotlib as mpl
import matplotlib.pyplot as plt


def simple1_interactive(browser, plotsections, plot_verbosity):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(np.linspace(1, 10, 100))

    plt.ioff()
    plt.show()
    return


if __name__ == "__main__":
    test_simple1(False, None, 10)
