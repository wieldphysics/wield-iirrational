#!/usr/bin/env python
# File: setup_example.py


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(
        'IIRrational',
        parent_package = parent_package,
        top_path       = top_path
    )
    #config.set_options(
    #    ignore_setup_xxx_py=True,
    #    assume_default_configuration=True,
    #    delegate_options_to_subpackages=True,
    #    quiet=True
    #)

    config.add_subpackage('lapack')
    return config



if __name__ == "__main__":
    from setuptools import setup
    from numpy.distutils.core import setup
    setup(
        name          = 'IIRrational',
        version       = '1.0.0.dev1',
        description   = "Fast Transfer Function Fitting",
        author        = "Lee McCuller",
        author_email  = "Lee.McCuller@gmail.com",
        py_modules    = ['IIRrational', 'IIRrational2'],
        configuration = configuration,
    )
