#!/usr/bin/env python
# -*- coding: utf-8 -*-
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: © 2021 Massachusetts Institute of Technology.
# SPDX-FileCopyrightText: © 2021 Lee McCuller <mcculler@caltech.edu>
# NOTICE: authors should document their contributions in concisely in NOTICE
# with details inline in source files, comments, and docstrings.
"""
"""
import argparse
import os
import sys
import numpy as np
import collections
from wield import declarative
from wield.bunch import Bunch, DeepBunch

try:
    from collections.abc import Mapping as MappingABC
except ImportError:
    from collections import Mapping as MappingABC

from .. import file_io

from .. import version as IIRversion
#from .. import auto_version as IIRversion_auto
from ..utilities.strings import padding_remove
from .. import fitters_ZPK

from . import arguments
from .arguments import pyargparse
from .arguments import cmdline
from . import fit_aid

from .data2filter import data2filter

# from .arguments import ArgumentError


class DataError(ValueError):
    pass


def find_groups(fdict, groups):
    gdict = dict()
    used_groups = []
    for key in reversed(groups):
        try:
            subgroup = file_io.subkey_search(fdict, key, default=declarative.NOARG)
        except TypeError:
            continue
        if subgroup is not declarative.NOARG:
            used_groups.append(used_groups)
            gdict.update(subgroup)
    return gdict, used_groups


def dump_obj(d):
    if isinstance(d, np.ndarray):
        if d.shape == ():
            return "[]"
        if len(d) < 8:
            return "[" + (", ".join(dump_obj(o) for o in d)) + "]"
        else:
            return "[{} {} array]".format(d.shape, d.dtype)
    elif isinstance(d, list):
        d = np.asarray(d)
        if len(d) < 8:
            return "[" + (", ".join(dump_obj(o) for o in d)) + "]"
        else:
            return "[{} {} list]".format(len(d), d.dtype)
    elif isinstance(d, tuple):
        d = np.asarray(d)
        if len(d) < 8:
            return "(" + (", ".join(dump_obj(o) for o in d)) + ")"
        else:
            return "[{} {} tuple]".format(len(d), d.dtype)
    else:
        return str(d)


def dump_fdict_keys(d, k=None, depth=-1, file=sys.stdout):
    subprint = False
    if isinstance(d, MappingABC):
        if k is not None:
            print("{}{}:".format("  " * depth, k), file=file)
        for subk, v in sorted(d.items()):
            dump_fdict_keys(v, subk, depth=depth + 1)
    elif isinstance(d, (np.ndarray, list, tuple)):
        d = np.asarray(d)
        if d.shape == ():
            return dump_fdict_keys(d.item(), k, depth, file=file)
        else:
            print("{}{}: {}".format("  " * depth, k, dump_obj(d)), file=file)
        if d.dtype == object and len(d) < 8:
            subprint = True
    else:
        print("{}{}: {}".format("  " * depth, k, str(d)), file=file)

    if subprint:
        for idx, v in enumerate(d):
            dump_fdict_keys(v, idx, depth=depth + 1)


def find_elements(fdict, groups, allow_dicts):
    used_keys = []
    elements = []
    for key in reversed(groups):
        try:
            val = file_io.subkey_search(fdict, key)
        except TypeError:
            continue
        except KeyError:
            continue

        if not allow_dicts and isinstance(val, MappingABC):
            continue

        used_keys.append(used_keys)
        elements.append(val)
    return elements, used_keys


def config_load(aid, args, kwfull, mode):
    # CONFIG File Loading
    cfiles = arguments.grab_kwargs(aid, args, cmdline.kw_hints, "config")
    cgroup = arguments.grab_kwargs(aid, args, cmdline.kw_hints, "config_group")

    fdicts_conf = dict()
    for cfname in reversed(cfiles):
        try:
            fB_conf = file_io.determine_type(cfname)
        except KeyError:
            raise arguments.ArgumentError(
                ("Unrecognized file type for {}").format(cfname)
            )
        if fB_conf.ftype == "csv":
            raise arguments.ArgumentError("Configuration files may not be csv")
        elif fB_conf.ftype == "ini":
            raise arguments.ArgumentError("ini files not supported")
        cgroup_use = cgroup
        if fB_conf.subkey is not None:
            cgroup_use = cmdline.mapcheck_cslist(fB_conf.subkey)

        fdict_conf = file_io.load_any(
            fname=fB_conf.fname,
            ftype=fB_conf.ftype,
        )
        fdicts_conf[cfname] = fdict_conf
        # detect and warn if more than one are found
        conf_dict, used_cgroups = find_groups(fdict_conf, cgroup_use)
        if len(used_cgroups) == 0:
            print(
                ("No config groups found in config file {}").format(fB_conf.fname),
                file=sys.stderr,
            )
            print("The config file layout is:", file=sys.stderr)
            dump_fdict_keys(fdict_conf, file=sys.stderr)
            # TODO dump data file layout
            raise SystemExit(1)
        elif len(used_cgroups) > 1:
            print(
                (
                    "Warning: multiple possible config groups found in config file {}"
                ).format(fB_conf.fname),
                file=sys.stderr,
            )
            print("Groups aggregated: {}".format(used_cgroups), file=sys.stderr)

        if mode == "printconfig":
            print("config file layout for mode=printconfig")
            print("file: {}".fortmat(fB_conf.fname))
            dump_fdict_keys(fdict_conf)

        kwfull.update(conf_dict)

    kwfull.update(args)
    return Bunch(
        cfiles=cfiles,
        fdicts=fdicts_conf,
        cgroup=cgroup,
    )


def data_load(aid, args, kwfull, mode, confB):
    # DATA FILE LOADING
    dfile = args["datafile"]
    try:
        fB_data = file_io.determine_type(dfile)
    except KeyError:
        raise arguments.ArgumentError(
            ("Unrecognized file type for data file {}").format(dfile)
        )

    dgroup = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "data_group")

    conf_dict = None
    if fB_data.ftype == "ini":
        raise arguments.ArgumentError("ini type not supported for data")
    if fB_data.ftype == "csv":

        if fB_data.subkey is None:
            parse_str = ",F,Xr,Xi,S,E"
        else:
            parse_str = fB_data.subkey

        F_list = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "frequency")
        Xr_list = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "dataXreal")
        Xi_list = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "dataXimaginary")
        SNR_list = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "dataSNR")
        COH_list = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "dataCOH")
        E_list = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "dataEmphasis")
        parse_map = {
            "F": F_list[0],
            "Xr": Xr_list[0],
            "Xi": Xi_list[0],
            "S": SNR_list[0],
            "C": COH_list[0],
            "E": E_list[0],
            ".": None,
            "*": None,
            "?": None,
            "": None,
        }

        try:
            fdict_data = file_io.load_csv(
                fname=fB_data.fname,
                parse_str=parse_str,
                parse_map=parse_map,
            )
        except KeyError as e:
            k = e.args[0]
            raise arguments.ArgumentError(
                (
                    "CSV parse string '{}' has unrecognized key {}. Keys must be within"
                    " ['F', 'Xr', 'Xi', 'S', 'E', '.']"
                ).format(parse_str, k)
            )
        data_dict = fdict_data

        if mode == "printdata":
            dump_fdict_keys(data_dict, file=sys.stderr)
            raise SystemExit(0)
    elif fB_data.ftype == "special":
        fdict_data = file_io.load_any(
            fname=fB_data.fname,
            ftype=fB_data.ftype,
        )
        data_dict = fdict_data
        conf_dict = fdict_data
    else:
        fdict_data = file_io.load_any(
            fname=fB_data.fname,
            ftype=fB_data.ftype,
        )

        if mode == "printdata":
            print("data file layout for mode=printdata")
            dump_fdict_keys(fdict_data)
            raise SystemExit(0)

        dgroup_use = dgroup
        if fB_data.subkey is not None:
            dgroup_use = cmdline.mapcheck_cslist(fB_data.subkey)

        if len(dgroup_use) == 0:
            data_dict = fdict_data
            used_dgroups = [None]
        else:
            # detect and warn if more than one are found
            data_dict, used_dgroups = find_groups(fdict_data, dgroup_use)

        if len(used_dgroups) == 0:
            print(
                ("No data groups found in data file {}").format(fB_data.fname),
                file=sys.stderr,
            )
            print("The data file layout is:", file=sys.stderr)
            dump_fdict_keys(fdict_data, file=sys.stderr)
            # TODO dump data file layout
            raise SystemExit(1)
        elif len(used_dgroups) > 1:
            print(
                ("Warning: multiple possible data groups found in data file {}").format(
                    fB_data.fname
                ),
                file=sys.stderr,
            )
            print("Groups aggregated: {}".format(used_dgroups), file=sys.stderr)

        # detect and warn if more than one are found
        conf_dict, used_cgroups = find_groups(fdict_data, confB.cgroup)
        if len(used_cgroups) > 1:
            print(
                (
                    "Warning: multiple possible config groups found in data file {}"
                ).format(fB_data.fname),
                file=sys.stderr,
            )
            print("Groups aggregated: {}".format(used_cgroups), file=sys.stderr)

        if mode == "printconfig":
            print("data file layout for mode=printconfig, ignoring data groups")
            print("file: {}".fortmat(fB_data.fname))
            # pop out the used data groups. Only OK since we are exiting after
            # the print
            for dg in used_dgroups:
                fdict_data.pop(dg, None)
            dump_fdict_keys(fdict_data)
            raise SystemExit(0)

    ##debugging
    # dump_fdict_keys(data_dict, file = sys.stderr)

    # if any were loaded and conf_dict is not None:
    if conf_dict:
        kwfull.update(conf_dict)
        # args always takes precedent
        kwfull.update(args)

    return Bunch(
        fdict=fdict_data,
        data_dict=data_dict,
    )


def data_assemble(aid, kwfull, dataB):
    def pull_data_element(elname):
        """
        pulls a single element from a list of possible elements
        """
        clist = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, elname)
        val, used = find_elements(dataB.data_dict, clist, allow_dicts=False)
        if len(used) > 1:
            raise DataError(
                ("multiple possible redundant elements in data group: {}").format(used)
            )
        elif len(used) > 0:
            return val[0], used[0], clist
        else:
            return None, None, clist

    def pull_data_elements(elname):
        """
        pulls a multiple element from a list of possible elements
        """
        clist = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, elname)
        val, used = find_elements(dataB.data_dict, clist, allow_dicts=False)
        if len(used) > 0:
            return val, used, clist
        else:
            return None, None, clist

    xfer = None
    try:
        XF, XF_used, XF_list = pull_data_element("frequency")
        if XF is None:
            raise DataError(
                (
                    "Frequency must be specified as one of {} (set by --frequency)"
                ).format(XF_list)
            )

        Xc, Xc_used, Xc_list = pull_data_element("dataXcomplex")
        if Xc is not None:
            xfer = np.array(Xc)
            if not np.iscomplexobj(xfer):
                raise (
                    RuntimeError(
                        "Xfer found, but does not appear to be complex {}".format(xfer)
                    )
                )

        Xr, Xr_used, Xr_list = pull_data_element("dataXreal")
        Xi, Xi_used, Xi_list = pull_data_element("dataXimaginary")
        if Xr is not None and Xi is not None:
            if xfer is not None:
                raise DataError("redundant elements in data group")
            xfer = np.asarray(Xr) + 1j * np.asarray(Xi)
        elif Xr is not None:
            raise DataError("Data {} given but not imag part".format(Xr_used))
        elif Xi is not None:
            raise DataError("Data {} given but not real part".format(Xi_used))

        Xa, Xa_used, Xa_list = pull_data_element("dataXamplitude")
        Xdb, Xdb_used, Xdb_list = pull_data_element("dataXdB")
        if Xa is not None and Xdb is not None:
            raise DataError("redundant amplitude elements in data group")
        elif Xdb is not None:
            Xa = 10 ** (np.asarray(Xdb) / 20.0)

        Xp, Xp_used, Xp_list = pull_data_element("dataXphaseRad")
        Xd, Xd_used, Xd_list = pull_data_element("dataXphaseDeg")
        if Xp is not None and Xd is not None:
            raise DataError("redundant phase elements in data group")
        elif Xd is not None:
            Xp = np.asarray(Xd) * np.pi / 180

        if Xa is not None and Xp is not None:
            if xfer is not None:
                raise DataError("redundant elements in data group")
            xfer = np.asarray(Xa) * np.exp(1j * np.asarray(Xp))
        elif Xa is not None:
            raise DataError("Data amplitude given but not phase part")
        elif Xp is not None:
            raise DataError("Data phase given but not amplitude part")

        if xfer is None:
            raise DataError("No transfer function elements given in data group")

        XS, XS_used, XS_list = pull_data_elements("dataSNR")
        if XS is not None:
            if len(XS) > 1:
                XS = tuple(np.asarray(_) for _ in XS)
            else:
                XS = XS[0]
        XC, XC_used, XC_list = pull_data_elements("dataCOH")
        if XC is not None:
            if len(XC) > 1:
                XC = tuple(np.asarray(_) for _ in XC)
            else:
                XC = XC[0]
        XE, XE_used, XE_list = pull_data_elements("dataEmphasis")
        if XE is not None:
            if len(XE) > 1:
                XE = tuple(np.asarray(_) for _ in XE)
            else:
                XE = XE[0]
    except DataError as e:
        print("Error: {}".format(e.args[0]))
        print("Data group layout is: ")
        dump_fdict_keys(dataB.data_dict)
        print("Error: {}".format(e.args[0]))
        raise SystemExit(1)

    data_dict = dict(
        xfer=xfer,
        F_Hz=XF,
        emphasis=XE,
        SNR=XS,
        COH=XC,
    )
    return data_dict


def check_pzargs(args):
    p = args.get("poles", None)
    z = args.get("zeros", None)
    P = args.get("poles_overlay", None)
    Z = args.get("zeros_overlay", None)
    msg = (
        "Error, must specify some poles if -{0} flag given, "
        " could be parse ambiguity with '-', use -{0}=-#.##,-#.##... "
        " or -{0} -#.## -#.## -#.## "
    )
    if p is not None and len(p) == 0:
        print(msg.format("p"), file=sys.stderr)
        SystemExit(1)
    elif z is not None and len(z) == 0:
        print(msg.format("z"), file=sys.stderr)
        SystemExit(1)
    elif P is not None and len(P) == 0:
        print(msg.format("P"), file=sys.stderr)
        SystemExit(1)
    elif Z is not None and len(Z) == 0:
        print(msg.format("Z"), file=sys.stderr)
        SystemExit(1)


def IIRrationalV2fit(argv=None, **kwargs):
    """ """
    if argv is None:
        argv = sys.argv[1:]
    elif isinstance(argv, str):
        import shlex

        argv = shlex.split(argv)


    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=False,
    )
    ap.add_argument('--cmdline_file', dest='cmdline_file', default=None)
    args, unknown_args = ap.parse_known_args(argv)
    if args.cmdline_file is not None:
        # print("COMMAND_LINE FILE: ", args.cmdline_file)
        # print(unknown_args)
        argv = unknown_args
        try:
            with open(args.cmdline_file) as F:
                cmdline_add = F.read()
                print("FILE: ", cmdline_add)
                import shlex
                argv2 = shlex.split(cmdline_add)
                print("argv2", argv2)
                argv = argv + argv2
        except FileNotFoundError:
            print("WARNING: file not found ", args.cmdline_file)
            pass

    aid = fit_aid.FitAid()

    kw_hints = dict(arguments.kw_hints)
    for name in arguments.standardargs.kw_hints_data.keys():
        del kw_hints[name]
    kw_hints.update(cmdline.kw_hints)

    ap = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
    )
    ap = pyargparse.kwdict_argparse(ap, kw_hints, groups_kw=cmdline.groups_kw)
    args = ap.parse_args(argv)
    args = vars(args)

    mode = arguments.grab_kwargs(aid, args, cmdline.kw_hints, "mode")
    kwfull = dict()

    confB = config_load(
        aid=aid,
        kwfull=kwfull,
        args=args,
        mode=mode,
    )

    ofile = args["outputfile"]
    if ofile == "":
        ofile = None
    overwrite = arguments.grab_kwargs(aid, args, cmdline.kw_hints, "overwrite")
    if ofile is not None:
        try:
            fB_output = file_io.determine_type(ofile)
        except KeyError:
            raise arguments.ArgumentError(
                ("Unrecognized file type for output file {}").format(ofile)
            )
        if not overwrite and os.path.exists(fB_output.fname):
            print(
                ("Error: output file {} already exists").format(ofile), file=sys.stderr
            )
            raise SystemExit(1)

    dataB = data_load(
        aid=aid,
        kwfull=kwfull,
        args=args,
        mode=mode,
        confB=confB,
    )

    if mode == "printsettings":
        print("settings dictionary after all loaded configurations")
        dump_fdict_keys(kwfull)
        raise SystemExit(0)

    # pull these early to make sure they are there and parsed before doing more work
    choose = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "choose")
    information = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "information")
    LIGO_foton = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "LIGO_foton")
    refine = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "refine")
    refine_file = arguments.grab_kwargs(aid, kwfull, cmdline.kw_hints, "refine_file")
    plot_fit = arguments.grab_kwargs(aid, args, cmdline.kw_hints, "plot_fit")
    plot_order = arguments.grab_kwargs(aid, args, cmdline.kw_hints, "plot_order")

    data_dict = data_assemble(
        aid=aid,
        kwfull=kwfull,
        dataB=dataB,
    )

    # Check these for error handing, since these can be zero for bad reasons
    check_pzargs(args=args)

    kwuse = dict(
        mode=mode,
    )
    arguments.transfer_kw(kwfull, kwuse, arguments.kw_hints)
    kwconfig = dict(kwuse)

    kwuse.update(data_dict)

    # TODO verbosity?
    print("Arguments: ")
    dump_fdict_keys(kwuse, depth=0)

    results = data2filter(
        # coding_map=fitters_ZPK.codings_s.coding_maps.SOS,
        # coding_map=fitters_ZPK.codings_s.coding_maps.SOSnl,
        **kwuse
    )

    version_dict = dict(
        data2filter="v2",
        version=IIRversion,
        # git_branch=IIRversion_auto.git_branch,
        # git_shorthash=IIRversion_auto.git_shorthash,
        # git_longhash=IIRversion_auto.git_longhash,
        # git_last_tag=IIRversion_auto.git_last_tag,
    )

    if plot_order is not None:
        try:
            print("plotting for -R, --plot_order={}".format(plot_order))
            results.investigate_order_plot(fname=plot_order)
        except Exception as E:
            print("Error: plotting exception {}".format(E), file=sys.stderr)

    if choose in ["shell", "prompt"]:
        choice_shell(
            results=results,
            data=DeepBunch(data_dict),
            config=DeepBunch(kwconfig),
            dfile=DeepBunch(dataB.fdict),
            cfiles=DeepBunch(confB.fdicts),
            version=DeepBunch(version_dict),
        )
    elif choose.startswith("auto"):
        remainder = choose[4:]
        if not remainder:
            order = 10
        else:
            try:
                order = int(remainder)
            except ValueError:
                print(
                    ("Error: auto order choice not an int, using auto10"),
                    file=sys.stderr,
                )
                order = 10
        order_arr = results.investigate_order_arrays().order
        num_lower = np.count_nonzero(order_arr < results.order)
        if num_lower > 0 and results.order > order:
            print(
                (
                    "Baseline order is {} > --choose=auto{}, and {} fit(s) are lower order, using shell"
                ).format(results.order, order, num_lower)
            )
            choice_shell(
                results=results,
                data=DeepBunch(data_dict),
                config=DeepBunch(kwconfig),
                dfile=DeepBunch(dataB.fdict),
                cfiles=DeepBunch(confB.fdicts),
                version=DeepBunch(version_dict),
            )
        else:
            print(
                (
                    "Baseline order is {} < --choose=auto{}, and {} fit(s) are lower order, using baseline"
                ).format(results.order, order, num_lower)
            )
    elif choose == 'best':
        results.choose(1000)
    elif choose.startswith("baseline"):
        remainder = choose[8:]
        if not remainder:
            order = None
        else:
            try:
                order = int(remainder)
            except ValueError:
                print(
                    ("Error: baseline order choice not an int, using baseline10"),
                    file=sys.stderr,
                )
                order = 10
        print("CHOOSING:", order)
        if order is not None and results.order > order:
            results.choose(order)
    else:
        try:
            order = int(choose)
        except ValueError:
            print(
                "Error: order choice not an int, using baseline order", file=sys.stderr
            )
        else:
            results.choose(order)

    if LIGO_foton in ["-", "None", "none"]:
        pass
    elif LIGO_foton == "Sf":
        print("-L, --LIGO_foton=Sf output:")
        print("%%%")
        print(results.as_foton_str_ZPKsf())
        print("%%%")
    else:
        print("Error: Unrecognized foton output type", file=sys.stderr)

    if refine:
        print(
            "%%% add or adjust the following poles and zeros to bootstrap the fit and refine. Disable this message with --no-refine"
        )
        print(results.as_refine_str())
        print("%%%")

    if refine_file:
        with open(refine_file, 'w') as F:
            F.write(results.as_refine_str())

    if plot_fit is not None:
        try:
            print("plotting for -r, --plot_fit={}".format(plot_fit))
            results.investigate_fit_plot(fname=plot_fit)
        except Exception as E:
            print("Error: plotting exception {}".format(E), file=sys.stderr)
            raise

    # TODO, need to ensure that this output form meets the common exchange
    # requirements for writing to the many filetypes
    output = collections.OrderedDict()
    output["fit"] = results.choice_export_dict()
    output["config"] = kwconfig
    output["fit_alternatives"] = results.alternatives_export_dict()
    output["data"] = data_dict
    output["version"] = version_dict
    # output['config_full']      = results.fit_aid.hints
    if information is not None:
        output["information"] = information

    if ofile is not None:
        file_io.write_any(
            fname=fB_output.fname,
            ftype=fB_output.ftype,
            fdict=output,
        )
    return results


def choice_shell(**kw):
    # First import the embed function
    from IPython.terminal.embed import InteractiveShellEmbed

    banner = padding_remove(
        """
    ------------------------------------------
    wield.iirrational fit choose investigation shell
    ------------------------------------------
    defined variables:
        results: result object of the wield.iirrational data2filter fit
        data:    loaded data dictionary
        config:  loaded config dictionary
        dfile:   original data file, normalized to dictionary form
        cfiles:  original config files, normalized to dictionary form
        exit:    function call exit() to leave shell

    use -j, --choose option to prevent dropping into this shell

    tab complete at ">results." to get methods of the results object
    to inspect the fit.

    typical usage:
    >results.investigate_order_plot()
    >results.choose(6)
    >exit()
    """
    )

    user_ns = dict(kw)

    ipshell = InteractiveShellEmbed(
        banner1=banner,
        user_ns=user_ns,
    )
    ipshell.enable_matplotlib()
    ipshell()


main = IIRrationalV2fit

if __name__ == "__main__":
    import sys

    IIRrationalV2fit()
