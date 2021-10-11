# -*- coding: utf-8 -*-
"""
"""
from __future__ import division, print_function, unicode_literals
import os.path as path
import subprocess


def pandoc_pdf(out_name):
    def exporter(md_name, folder):
        subprocess.call(['pandoc', '-o', out_name, md_name], cwd = folder)
    return exporter

def npm_markdown_pdf(out_name):
    def exporter(md_name, folder):
        subprocess.call(['markdown-pdf', md_name, '-o', out_name], cwd = folder)
    return exporter


