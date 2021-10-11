# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import os
import re
import errno
import os.path as path

from ..utilities.mpl import asavefig
from ..utilities.strings import padding_remove

def header_line(
        level = 0,
        section_history = (),
        heading = ''
):
    if not heading:
        heading = ''
    section_link = '.'.join([str(n) for n in section_history])
    if section_link:
        section_link = ' `{0}` '.format(section_link)
    else:
        section_link = ' '
    return ('#' * (level + 1)) + section_link + heading

def docdb_md_section(
    F_md,
    doc_db,
    section_history,
    heading_level,
    md_verbosity    = None,
    plot_verbosity  = None,
    plot_fname_func = None,
    do_plot         = True,
    regex_secname   = '',
    use_ext         = True,
    plot_format     = 'png',
):
    if md_verbosity is not None:
        if doc_db.verbosity:
            if doc_db.verbosity > md_verbosity:
                return

    header = header_line(
        heading_level,
        section_history,
        doc_db.name,
    )

    lines = [
        header,
        '',
    ]

    if doc_db.method:
        lines.append('> `{0}`'.format(padding_remove(doc_db.method)))
        lines.append('')

    if doc_db.about:
        lines.append(padding_remove(doc_db.about))
        lines.append('')

    if not do_plot:
        if re.findall(regex_secname, doc_db.name):
            do_plot = True

    will_plot = True
    if plot_verbosity is not None:
        if doc_db.verbosity:
            if doc_db.verbosity > plot_verbosity:
                will_plot = False

    if do_plot and will_plot and doc_db.plotter and doc_db.fitter:
        title = doc_db.plot_title
        if not title:
            title = ''
        plot_fname = plot_fname_func(
            section_history,
            doc_db,
            with_folder = False,
        )
        if use_ext:
            ext = '.' + plot_format
        else:
            ext = ''
        lines.append('![{title}]({0})'.format(plot_fname + ext, title = title))
        lines.append('')

    F_md.writelines((l + '\n' for l in lines))

    if doc_db.sequence:
        for idx, seq_name in enumerate(doc_db.sequence):
            docdb_md_section(
                F_md,
                doc_db[seq_name],
                section_history = section_history + [idx + 1],
                heading_level   = heading_level + 1,
                md_verbosity    = md_verbosity,
                plot_verbosity  = plot_verbosity,
                plot_fname_func = plot_fname_func,
                do_plot         = do_plot,
                regex_secname   = regex_secname,
                plot_format     = plot_format,
                use_ext         = use_ext,
            )
    return

def docdb_md_plots(
    doc_db,
    section_history,
    plot_verbosity  = None,
    plot_fname_func = None,
    do_plot         = True,
    regex_secname   = '',
    clear_plots     = False,
    plot_format     = 'png',
    dpi             = 100,
):
    if plot_verbosity is not None:
        if doc_db.verbosity:
            if doc_db.verbosity > plot_verbosity:
                return

    if not do_plot:
        if re.findall(regex_secname, doc_db.name):
            do_plot = True

    if do_plot and doc_db.plotter and doc_db.fitter:
        plot_fname = plot_fname_func(
            section_history,
            doc_db,
            with_folder = True,
        )
        print("PLOTTING: ", plot_fname)
        Fb = doc_db.plotter(doc_db.fitter)
        Fb.formats[plot_format].use = True
        Fb.formats[plot_format].facecolorize = True
        Fb.formats[plot_format].dpi = dpi
        Fb.save_show = False
        Fb.save(plot_fname, fixname = False)
    elif clear_plots:
        plot_fname = plot_fname_func(
            section_history,
            doc_db,
            with_folder = True,
        )
        print("REMOVING: ", plot_fname)
        try:
            os.remove(plot_fname + '.' + plot_format)
        except OSError as e:
            if e.errno != errno.ENOENT:
                print(e)

    if doc_db.sequence:
        for idx, seq_name in enumerate(doc_db.sequence):
            docdb_md_plots(
                doc_db[seq_name],
                section_history = section_history + [idx + 1],
                plot_verbosity  = plot_verbosity,
                plot_fname_func = plot_fname_func,
                do_plot         = do_plot,
                regex_secname   = regex_secname,
                clear_plots     = clear_plots,
                plot_format     = plot_format,
                dpi             = dpi,
            )
    return

def docdb2markdown(
    folder_name,
    doc_db_root,
    md_name            = 'digest.md',
    include_plots      = True,
    include_md         = True,
    plot_verbosity     = None,
    md_verbosity       = None,
    clear_plots        = False,
    regex_plotsections = None,
    ipy_display        = False,
    use_ext            = True,
    plot_format        = 'png',
    dpi                = 100,
    exporter           = None,
    MP_workers         = 1,
):
    try:
        os.makedirs(folder_name)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    def plot_fname_func(
            section_history,
            doc_db,
            with_folder = True,
    ):
        section_link = '-'.join([str(n) for n in section_history])
        if section_link:
            fname = 'plot-' + section_link
        else:
            fname = 'plot-main'
        fname = fname.replace('_', '-')
        if with_folder:
            return path.join(folder_name, fname)
        else:
            return fname

    if regex_plotsections is not None:
        do_plot = False
    else:
        do_plot = True

    if include_md:
        with open(path.join(folder_name, md_name), 'w') as F_md:
            docdb_md_section(
                F_md,
                doc_db_root,
                section_history = [],
                heading_level   = 0,
                md_verbosity    = None,
                plot_fname_func = plot_fname_func,
                plot_verbosity  = plot_verbosity,
                regex_secname   = regex_plotsections,
                do_plot         = do_plot,
                plot_format     = plot_format,
                use_ext         = use_ext,
            )

    if include_plots:
        with asavefig.pool(workers = MP_workers):
            docdb_md_plots(
                doc_db_root,
                section_history = [],
                plot_verbosity  = plot_verbosity,
                plot_fname_func = plot_fname_func,
                regex_secname   = regex_plotsections,
                do_plot         = do_plot,
                clear_plots     = clear_plots,
                plot_format     = plot_format,
                dpi             = dpi,
            )

    if ipy_display:
        with open(path.join(folder_name, md_name), 'r') as F_md:
            md_data = F_md.read()
        md_data = re.sub(r'\((.*)\.png\)', r'({0}/\1.png)'.format(folder_name), md_data)
        try:
            from IPython.display import display_markdown
        except ImportError:
            print("Can't Import IPython.display")
        else:
            display_markdown(md_data, raw = True)

    if exporter is not None:
        exporter(
            md_name = md_name,
            folder = folder_name,
        )
    return
