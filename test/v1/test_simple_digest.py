# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import pytest

import os.path as path

from IIRrational.testing import IIRrational_data
from IIRrational import v1
from IIRrational.annotate.exporters import npm_markdown_pdf


def test_simple_digest(browser, plotsections, plot_verbosity):
    data_name = 'simple2'
    out = v1.data2filter(
        IIRrational_data(
            data_name,
            instance_num = 1,
            set_num = 2,
        ),
        order_initial = 30,
        delay_s = 0,
    )

    folder = path.join(path.split(__file__)[0], 'out-' + data_name)
    print('FOLDER', folder)

    if browser:
        out.digest_write(
            folder,
            md_name            = 'browser.md',
            plot_verbosity     = 10,
            clear_plots        = True,
            regex_plotsections = plotsections,
            plot_format        = 'png',
            dpi                = 100,
        )
    if plot_verbosity is None:
        plot_verbosity = 4
    out.digest_write(
        folder,
        plot_verbosity     = plot_verbosity,
        clear_plots        = True,
        regex_plotsections = plotsections,
        plot_format        = 'png',
        dpi                = 200,
        exporter           = npm_markdown_pdf('digest.pdf'),
        MP_workers         = 1,
    )
    return


if __name__=='__main__':
    test_simple1(False, None, 10)
