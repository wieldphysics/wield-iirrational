"""
(note, this is generated from the command line version, 80% representative
of the options for the direct function call)
optional arguments:
  -h, --help            show this help message and exit

  data.file             Specify data file. This file may be .h5, .hdf, .mat, .csv, .txt.
                        Can specify or override the filetype and internal format using a
                        namespace syntax '::'.
                         - "xfer.___::hdf" for hdf5 with an unrecognized suffix
                         - "xfer.out::csv::,F,Xr,Xi,W" to interpret as csv with 4 columns including
                           data real and imaginary parts
                         - "xfer.out::csv::-F,Xm,Xp,W" to interpret as csv with 4 columns whitespace
                           separated including data magnitude and phase [radians]
                        all of the CSV column specifiers are:
                         - first specifier is the delimiter, one of any puncuation such as  ",:|".
                           the symbol '-' indicates to use whitespace of any kind.
                         - F, f for frequency in [Hz]
                         - Xr, Xi for xfer real and imaginary parts
                         - Xm, Xa for xfer magnitude or equivalently amplitude
                         - Xdb, XdB xfer in dBs
                         - Xp for xfer phase in radians
                         - Xd for xfer phase in degrees
                         - S, W, w for weights or SNR
                         - E, e for emphasis
                         - . for ignore column
                        by default, the CSV is expected to be formatted as ",XrXiFSE". If columns
                        are missing they are dropped
                        matlab and h5 structures are assumed to have keys:
                         'F_Hz', 'xfer', 'SNR', 'emphasis',
                        but arguments in the "data" group can override these keys, or specify
                        columns.
  output.file           Specify the output file to store fit results including zpk of the
                        chosen and alternative fits, the configurations, and versioning information.
                        Possible output extensions are .h5, .hdf, .mat, .pkl, .json, .yaml. Binary
                        formats .mat, .h5, .pkl will include the original data, for full reconstruction
                        of the fit. IIRrational may be called on the output file to rerun the fit.
  -c CONFIG, --config CONFIG
                        Specify configuration file. This file may be .json, .yaml, .ini, .h5, .mat.
                        These will be interpreted as a dictionary or structure full of keyword
                        arguments. This may be specified multiple times to overlay settings.
                        Additional command line setting always override these configuration files.
  -L {-,none,Sf}, --LIGO_foton {-,none,Sf}
                        What type of foton output would you prefer?
                          - '-' or 'None' to suppress foton output
                          - Sf for frequency domain (default)
                          - n for normalized (not yet supported)
                          - Sw for radial (not yet supported)
  -j CHOOSE, --choose CHOOSE
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
  -m {full,fit,reduce,rational,dumpargs,dumpargs_full,printdata,printconfig,printsettings}, --mode {full,fit,reduce,rational,dumpargs,dumpargs_full,printdata,printconfig,printsettings}
                        Fitting mode to use, to change the automation level. Must be one of
                         - "full": to use rational fitting to get initial parameters, then alternate
                                   optimization and order reduction
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
  -N F_NYQUIST_HZ, --F_nyquist_Hz F_NYQUIST_HZ
                        This selects the Nyquist frequency for Z-domain fitting. If None or not
                        specified, the fits are done in the "Sf" domain - the s domain scaled for
                        units of frequency.
  -I INFORMATION, --information INFORMATION
                        A string to detail what this fit is of. For documentation purposes.
  --overwrite           Allow the output file to overwrite the previous output file (be careful).

groups:
  The data and config file formats (except csv) all load internally to a
  dictionary or structure representation. These arguments specify the
  keys or groups storing the relevant data or configuration subkeys. These
  also control the output groups when saving files.

  -C CONFIG_GROUP, --config_group CONFIG_GROUP
                        Group(s) within the loaded config or data file to search for
                        configuration elements. May be a comma separated list to aggregate groups.
                        defaults to "config,conf,IIRconfig,IIRconf". The first element will be
                        the group where configurations are stored to the output file.
  -D DATA_GROUP, --data_group DATA_GROUP
                        Group(s) within the loaded data file to search for data elements
                        (see "data" section). May be a comma separated list to aggregate groups.
                        defaults to "data,IIRdata,xfer". The first specified element will be
                        the group that data is stored to the output file.

data:
  Specify keys to search within the data group for specific arrays.
  -F,   --frequency,
  -Xc,  --dataXcomplex
  -Xr,  --dataXreal,      -Xi,  --dataXimaginary
  -Xa,  --dataXamplitude, -Xm,  --dataXmagnitude
  -XdB, --dataXdB,        -Xdb, --dataXdb
  -Xd,  --dataXphaseDeg,  -Xp,  --dataXphaseRad
  -S,   --dataSNR,        -W,   --dataW
  -E,   --dataEmphasis

order:
  --relative_degree RELATIVE_DEGREE
                        Sets the initial relative degree (number zeros minus polse).
                        Defaults to None, which will be the midpoint of the min and max if they
                        are specified, otherwise this will default to 0, which typically will
                        still fit filters to the correct degree. If constrained to be different
                        than the data, the filter will enter the asymptotic regime set by the
                        relative degree within 2x of the 'root_bandwidth_Hz_max' setting.
  --relative_degree_max RELATIVE_DEGREE_MAX
                        Maximum value for the filter relative degree (number zeros minus polse).
                        Defaults to the relative degree (which may be None). If this value is
                        None, the degree is unconstrained and the fit will land at some degree that fits well.
  --relative_degree_min RELATIVE_DEGREE_MIN
                        Minimum value for the filter relative degree (number zeros minus polse).
                        Defaults to the relative degree (which may be None). If this value is
                        None, the degree is unconstrained and the fit will land at some degree that fits well.
  --total_degree_min TOTAL_DEGREE_MIN
                        Minimum degree to search through during the successive order reduction phase.
                        Defaults to 2. If None, then successive reduction will not be performed.
  --order_initial ORDER_INITIAL
                        Order to use in rational fitting. If not specified, the order is
                        increased until the residuals plateau.
  --order_max ORDER_MAX
                        Maximium order that the rational fitter order estimation may use.
                        Increasing this value from the default may impact numerical stability.
  --order_min ORDER_MIN
                        Minimum order to use during order estimation. Smaller values will speed
                        up fits, but may fail to fix complex data.

plots:
  -r PLOT_FIT, --plot_fit PLOT_FIT
                        filename to plot (review) the chosen fit to. May be any extension supported by matplotlib.
                        [.pdf, .svg, .png, .jpg, ...]
  -R PLOT_ORDER, --plot_order PLOT_ORDER
                        filename to plot (Review) potential orders and residuals to. May be any extension supported by matplotlib.
                        [.pdf, .svg, .png, .jpg, ...]

delay:
  --delay_s_max DELAY_S_MAX
                        The maximum delay in seconds. Defaults to None, in which case delay is
                        never fit as a free parameter.
  --delay_s DELAY_S     Use this delay (in seconds) for the initial fitting up until the
                        "baseline" fit determination is complete. By default this is the same
                        as delay_s_min, unless delay_s_min is negative, in which case it defaults
                        to the smaller of 0 seconds or delay_s_max.
  --delay_s_min DELAY_S_MIN
                        The minimum delay in seconds. Defaults to 0. Also sets the typical
                        default value of the delay used for the initial part of the fits.
                        Must always be specified (cannot be None).

logging:
  -l LOG_LEVEL, --log_level LOG_LEVEL
                        Log level default for all logging types
  --log_level_alert LOG_LEVEL_ALERT
                        Data on useful statistical tests, particularly those which cause input
                        data/SNR/emphasis to be reinterpreted. Typically logs at level 2-4.
  --log_level_debug LOG_LEVEL_DEBUG
                        Debugging reports used for development.
  --log_level_info LOG_LEVEL_INFO
                        Logging on miscellaneous details.
  --log_level_rationale LOG_LEVEL_RATIONALE
                        Detailed explanations of tests and algorithms. Extremely verbose and
                        intended for first users and digest reports.
  --log_level_warn LOG_LEVEL_WARN
                        Warnings about failed tests, or data operating in a regime unexpected to
                        work or that validation should be made.

SNR Adjustments:
  --SNR_cuttoff SNR_CUTTOFF
                        Maximum SNR. Hard cutoff for SNR above this level.
  --SNR_estimate_width SNR_ESTIMATE_WIDTH
                        The window width of nearby points to use for averaging and median filtering
                        when estimating the SNR via windowed sample variance. Default is 10.
                        None or 0 indicates to not try.
  --SNR_regularize_ratio SNR_REGULARIZE_RATIO
                        Ratio of effective data points, determined by the dynamic range of the
                        SNR weighting, to the actual number of data points. This parameter
                        adjusts how the SNR is adjusted to ensure sufficient data points for a
                        good fit.
  --SNR_regularize_scale SNR_REGULARIZE_SCALE
                        Similar to SNR_regularize_ratio, but instead of a direct ratio, the ratio
                        is determined as 'SNR_regularize_scale' / max(SNR). For very high SNR
                        fits, the SNR should be whitened significantly using the ratio adjustment,
                        since the error is dominated by systematics rather than statistical noise.
                        The default is 10, causing a %90 coverage ratio when the highest SNR is
                        100 (1% statistical error).
  --trust_SNR           Overall parameter determining if SNR/W statistical weighting input should
                        be trusted. If False (the default) then additional statistical tests are
                        run to adjust the SNR for better fits. If the user can fully trust their
                        SNR estimators and fits to be unbiased, then the user likely does not
                        need this tool.

advanced:
  --resavg_RthreshOrdC RESAVG_RTHRESHORDC
                        The threshold relative change in the average residuals to accept a fit
                        of equal order. Only used for the "baseline" fit determination before
                        delay is activated and not used during the total order reduction.
  --resavg_RthreshOrdDn RESAVG_RTHRESHORDDN
                        The threshold relative change in the average residuals to accept a fit
                        of lower order. Only used for the "baseline" fit determination before
                        delay is activated and not used during the total order reduction.
  --resavg_RthreshOrdUp RESAVG_RTHRESHORDUP
                        The threshold relative change in the average residuals to accept a fit
                        of higher order. Only used for the "baseline" fit determination before
                        delay is activated and not used during the total order reduction.

emphasis:
  --emphasis_order_initial EMPHASIS_ORDER_INITIAL
                        Order to use in rational fitting during the emphasis application stage.
  --emphasis_order_max EMPHASIS_ORDER_MAX
                        Maximium order that the rational fitter order estimation may use during
                        the emphasis application stage.
  --emphasis_order_min EMPHASIS_ORDER_MIN
                        Minimum order to use during order estimation during the emphasis application
                        stage.

fitting:
  -k GAIN, --gain GAIN  Initial gain for the fit
  -p [POLES [POLES ...]], --poles [POLES [POLES ...]]
                        Initial poles used in the fit.
  -P [POLES_OVERLAY [POLES_OVERLAY ...]], --poles_overlay [POLES_OVERLAY [POLES_OVERLAY ...]]
                        poles guaranteed to be part of the fit.
  -z [ZEROS [ZEROS ...]], --zeros [ZEROS [ZEROS ...]]
                        Initial zeros used in the fit.
  -Z [ZEROS_OVERLAY [ZEROS_OVERLAY ...]], --zeros_overlay [ZEROS_OVERLAY [ZEROS_OVERLAY ...]]
                        zeros guaranteed to be part of the fit.

operating mode:
  --suggest             How to use the rational fitter when provided an initial ZPK argument.
                        the default of False causes the rational to overlay the ZPK suggestion,
                        forcing the roots in the suggestion to be used during the optimization
                        stage. If set to True, the ZPK will suggest initial poles to use during
                        the fit, which can accellerate convergence or cause it to start at a
                        lower initial order, speeding up future fitting operations.

residuals:
  --alternate_residuals ALTERNATE_RESIDUALS
                        Flag to switch from using log/phase residuals to using dual
                        ratio residuals as the principle fit residuals. Log/phase is typically
                        better behaved.
  --residuals_type {log,dualA,dualB,poles,zeros}
                        Standard residuals type to use for the optimizations. Must be one of the
                        definitions below, where R=fit/data, the ratio of fit to data and W is
                        the SNR:
                            log:   W*(ln(|R|) + 1j*R.imag/|R|)
                            dualA: W*(R + 1/R - 2)/2
                            dualB: W*(R - 1/R)/2
                            poles: W*(1/R - 1)
                            zeros: W*(R - 1)
  --residuals_type_alt {log,dualA,dualB,poles,zeros}
                        Standard alternate residuals type to use during optimization annealing.
                        See residuals_type for options.

speed:
  --greedy_order GREEDY_ORDER
                        Does only partial optimization during the reduce step. Greatly speeds
                        up order reduction, but may not be as effective. Defaults to 30, where
                        it uses greedy optimization until reaching order 30, then uses combinatoric
                        optimization.

tunings:
  --distance_limit_scale DISTANCE_LIMIT_SCALE
                        Scaling for how to limit root bandwidth for roots in-between data points.
                        Increase to force roots to lower Q, decrease to allow higher Qs.
  --root_bandwidth_Hz_max ROOT_BANDWIDTH_HZ_MAX
                        Maximum bandwidth of any pole or zero in the fit. This affects the
                        asymptotic rolloff when relative_degree is constrained.
"""
#this file is auto generated during the setup.py build phase.

