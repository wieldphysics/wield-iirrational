#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@mit.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
from wavestate.bunch import Bunch

from .base import (
    ArgumentError,
    # mapcheck_bool,
    # cplx_iIjJ,
)

data_args_hide = True


def mapcheck_cslist(aid, aname, val):
    if isinstance(val, (list, tuple)):
        return val

    vals = val.split(",")
    vals = [v for v in [v.strip() for v in vals] if v != ""]
    return vals


def mapcheck_choose(aid, aname, val):
    """
    How to decide the optimal order. Special parameters are:
     - "prompt", "shell", "interactive" (the default) To enter a command prompt that plots and requests
       a choice
     - "auto10" the default choice to use the shell if the baseline order is
       above 10 unless there is no choice. (10 may be changed)
     - "baseline" to use the baseline fit before the ChiSq starts rising.
     - "baseline10" to use the baseline "best fit" unless it is over order 10
                    (10 can be configured to other integers)
     - integer. To force a choice of order. It will always use the best fit
       under this order.
    """
    modes = ["shell", "prompt", "interactive", "auto10", "baseline", "baseline10", "best"]
    try:
        if isinstance(val, str):
            val = val.lower()
            if val in modes:
                return val
            elif val.startswith("baseline"):
                remainder = val[8:]
                if int(remainder) < 0:
                    raise ValueError()
                return val
            
        # this runs the valueerror test that it can be parsed
        if int(val) < 0:
            raise ValueError()
        return val
    except ValueError:
        raise ArgumentError(
            (
                "Argument {}={} must be an +integer or one of {}".format(
                    aname, val, modes
                )
            )
        )


def mapcheck_print_(aid, aname, val):
    """
    How to decide the optimal order. Special parameters are:
     - "prompt", "shell", "interactive" (the default) To enter a command prompt that plots and requests
       a choice
     - "auto10" the default choice to use the shell if the baseline order is
       above 10 unless there is no choice. (10 may be changed)
     - "baseline" to use the baseline fit before the ChiSq starts rising.
     - "baseline10" to use the baseline "best fit" unless it is over order 10
                    (10 can be configured to other integers)
     - integer. To force a choice of order. It will always use the best fit
       under this order.
    """
    modes = ["shell", "prompt", "interactive", "auto10", "baseline", "baseline10"]
    try:
        if isinstance(val, str):
            val = val.lower()
            if val in modes:
                return val
        # this runs the valueerror test that it can be parsed
        if int(val) < 0:
            raise ValueError()
        return val
    except ValueError:
        raise ArgumentError(
            (
                "Argument {}={} must be an +integer or one of {}".format(
                    aname, val, modes
                )
            )
        )



kw_hints = Bunch(
    datafile=dict(
        APpositional=True,
        APpriority=-2,
        APmetavar="data.file",
        about="""
        Specify data file. This file may be .h5, .hdf, .mat, .yaml.
        these filetypes have internal directory structure and the data are by
        default assumed to be in the parameters:
         'F_Hz', 'xfer', 'SNR', 'emphasis', inside of the 'data' group.
        Arguments such is -Xc and -F and -D in the "data" options group can override
        these keys. It is possible also to use csv files, but the syntax is complicated.
        """,
    ),
    outputfile=dict(
        APpositional=True,
        APpriority=-2,
        APmetavar="output.file",
        about="""
        Specify the output file to store fit results including zpk of the
        chosen and alternative fits, the configurations, and versioning information.
        Possible output extensions are .h5, .hdf, .mat, .pkl, .json, .yaml. Binary
        formats .mat, .h5, .pkl will include the original data, for full reconstruction
        of the fit. wavestate.iirrational may be called on the output file to rerun the fit.
        """,
    ),
    overwrite=dict(
        APaction="store_true",
        default=False,
        about="""
        Allow the output file to overwrite the previous output file (be careful).
        """,
    ),
    mode=dict(
        APshort="-m",
        APpriority=3,
        APchoices=[
            "full",
            "full2x",
            "fullAAA",
            "onlyAAA",
            "AAA",
            "fit",
            "reduce",
            "rational",
            "rational2x",
            "dumpargs",
            "dumpargs_full",
            "printdata",
            "printconfig",
            "printsettings",
        ],
        default="full",
        about="""
        Fitting mode to use, to change the automation level. Must be one of
         - "full": to use rational fitting to get initial parameters, then alternate
                   optimization and order reduction
         - "full2x": runs the full optimization twice to refine it. Can often use a lower initial order.
         - "rational": to use rational fitting to get initial parameters and then
                       stop after the simplest order reduction. Useful to use with
                       other fitters. Does not perform delay fitting.
         - "fit": To take an initial ZPK guess, and only fit/optimize it. Good
                  for refining previous fits.
         - "reduce": To take an initial ZPK guess, and then alternate fitting and
                     order reduction.
         - "dumpargs": Dump arguments to determine current settings
         - "dumpargs_full": Dump arguments and all generated settings to see the full run setup
         - "printdata": Dump the layout of the data file specified
         - "printconfig": dump the layout of the config files specified
         - "printsettings": dump the fully loaded settings specified
        """,
    ),
    config=dict(
        APshort="-c",
        APaction="append",
        APpriority=0,
        default=[],
        about="""
        Specify configuration file. This file may be .json, .yaml, .ini, .h5, .mat.
        These will be interpreted as a dictionary or structure full of keyword
        arguments. This may be specified multiple times to overlay settings.
        Additional command line setting always override these configuration files.
        """,
    ),
    data_group=dict(
        APshort="-D",
        APgroup="groups",
        APpriority=1,
        default=["data", "IIRdata", "xfer"],
        mapcheck=mapcheck_cslist,
        about="""
        Group(s) within the loaded data file to search for data elements
        (see "data" section). May be a comma separated list to aggregate groups.
        defaults to "data,IIRdata,xfer". The first specified element will be
        the group that data is stored to the output file.
        """,
    ),
    config_group=dict(
        APgroup="groups",
        APpriority=1,
        default=["config", "conf", "IIRconfig", "IIRconf"],
        mapcheck=mapcheck_cslist,
        about="""
        Group(s) within the loaded config or data file to search for
        configuration elements. May be a comma separated list to aggregate groups.
        defaults to "config,conf,IIRconfig,IIRconf". The first element will be
        the group where configurations are stored to the output file.
        """,
    ),
    plot_fit=dict(
        APshort="-r",
        APgroup="plots",
        APpriority=10,
        default=None,
        about="""
        filename to plot (review) the chosen fit to. May be any extension supported by matplotlib.
        [.pdf, .svg, .png, .jpg, ...]
        """,
    ),
    plot_order=dict(
        APshort="-R",
        APgroup="plots",
        APpriority=10,
        default=None,
        about="""
        filename to plot (Review) potential orders and residuals to. May be any extension supported by matplotlib.
        [.pdf, .svg, .png, .jpg, ...]
        """,
    ),
    choose=dict(
        APshort="-j",
        APpriority=2,
        # TODO, type
        default="auto10",
        mapcheck=mapcheck_choose,
        about="""
        How to decide the optimal order. Special parameters are:
         - "prompt", "shell", "interactive" (the default) To enter a command prompt that plots and requests
           a choice
         - "auto10" the default choice to use the shell if the baseline order is
           above 10 unless there is no choice. (10 may be changed)
         - "baseline" to use the baseline fit before the ChiSq starts rising.
         - "baseline10" to use the baseline "best fit" unless it is over order 10
                        (10 can be configured to other integers)
         - integer. To force a choice of order. It will always use the best fit
           under this order.
        """,
    ),
    LIGO_foton=dict(
        APshort="-L",
        APpriority=2,
        APchoices=[
            "-",
            "none",
            "Sf",
        ],
        # TODO, type
        default="Sf",
        about="""
        What type of foton output would you prefer?
          - '-' or 'None' to suppress foton output
          - Sf for frequency domain (default)
          - n for normalized (not yet supported)
          - Sw for radial (not yet supported)
        """,
    ),
    refine=dict(
        APflags=["--no-refine"],
        APaction="store_false",
        default=True,
        APpriority=2,
        about="""
        disable the output hints of --zeros and --poles to run a second time and refine.
        """,
    ),
    refine_file=dict(
        APflags=["--refine_file"],
        default=None,
        APpriority=4,
        about="""
        indicate a file to write the hints of --zeros and --poles into.
        """,
    ),
    information=dict(
        APshort="-I",
        APpriority=6,
        default=None,
        about="""
        A string to detail what this fit is of. For documentation purposes.
        """,
    ),
    frequency=dict(
        APshort="-F",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "frequency_Hz",
            "freq_Hz",
            "F_Hz",
            "f_Hz",
            "frequency",
            "freq",
            "F",
            "f",
        ],
        about="key or column number to find complex data array",
    ),
    dataXcomplex=dict(
        APshort="-Xc",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "xfer",
            "data",
            "X",
            "transfer",
        ],
        about="key or column number to find complex data array",
    ),
    dataXreal=dict(
        APshort="-Xr",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "xfer.real",
            "xfer.r",
            "xfer_real",
            "xfer_r",
            "data.real",
            "data.r",
            "data_real",
            "data_r",
            "X.real",
            "X.r",
            "Xreal",
            "Xr",
            "transfer_real",
            "transfer.real",
            "transfer_r",
            "transfer.r",
        ],
        about="key or column number to find data array",
    ),
    dataXimaginary=dict(
        APshort="-Xi",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "xfer.imag",
            "xfer.i",
            "xfer_imag",
            "xfer_i",
            "X.imag",
            "X.i",
            "Ximag",
            "Xi",
            "transfer_imag",
            "transfer.imag",
            "transfer_i",
            "transfer.i",
        ],
        about="key or column number to find data array",
    ),
    dataXamplitude=dict(
        APflags=["-Xa", "-Xm", "--dataXamplitude", "--dataXmagnitude"],
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "xfer_mag",
            "xfer_magnitude",
            "xfer_amplitude",
            "xfer_amp" "Xa",
            "Xm",
            "Xamp",
            "Xmag" "amplitude",
            "magnitude",
            "amp",
            "mag",
        ],
        about="key or column number to find data array",
    ),
    dataXdB=dict(
        APflags=["-XdB", "-Xdb", "--dataXdB", "--dataXdb"],
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        default=["xfer_db", "xfer_dB", "Xdb", "XdB", "db", "dB"],
        about="key or column number to find data array",
    ),
    dataXphaseRad=dict(
        APshort="-Xp",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "xfer_phase",
            "xfer_radians",
            "xfer_rad",
            "Xrad",
            "Xphase",
            "phase",
            "radians",
            "rad",
        ],
        about="key or column number to find data array",
    ),
    dataXphaseDeg=dict(
        APshort="-Xd",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "xfer_degrees",
            "xfer_deg",
            "Xdeg",
            "Xdegrees",
            "degrees",
            "deg",
        ],
        about="key or column number to find data array",
    ),
    dataSNR=dict(
        APflags=["-W", "-S", "--dataSNR", "--dataWeights"],
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "SNR",
            "snr",
            "signal2noise" "weight",
            "Weight",
            "weights",
            "Weights",
            "W",
            "w",
        ],
        about="key or column number to find SNR array",
    ),
    dataCOH=dict(
        APflags=["-C", "--dataCOH"],
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        mapcheck=mapcheck_cslist,
        default=[
            "COH",
            "coh",
        ],
        about="""key or column numbers to find a coherence arrays. These are
        then used to compute the implied SNR. This argument can be specified
        in a comma-separated list to give multiple coherences to combine into
        a single SNR array.
        """,
    ),
    dataEmphasis=dict(
        APshort="-E",
        APhide=data_args_hide,
        APgroup="data",
        APpriority=5,
        APmetavar="key",
        default=["emphasis", "emph", "Emph", "E", "e"],
        mapcheck=mapcheck_cslist,
        about="key or column number to find data array",
    ),
)

groups_kw = dict(
    groups=dict(
        APpriority=4,
        about="""
        The data and config file formats (except csv) all load internally to a
        dictionary or structure representation. These arguments specify the
        keys or groups storing the relevant data or configuration subkeys. These
        also control the output groups when saving files.
    """,
    ),
    data=dict(
        APpriority=5,
        about="""
    Specify keys to search within the data group for specific arrays.
        -F,   --frequency,
        -Xc,  --dataXcomplex
        -Xr,  --dataXreal,      -Xi,  --dataXimaginary
        -Xa,  --dataXamplitude, -Xm,  --dataXmagnitude
        -XdB, --dataXdB,        -Xdb, --dataXdb
        -Xd,  --dataXphaseDeg,  -Xp,  --dataXphaseRad
        -S,   --dataSNR,        -W,   --dataW
        -E,   --dataEmphasis
    """,
    ),
)
