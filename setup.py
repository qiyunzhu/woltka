#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages

import woltka.__init__ as init


classes = """
    Development Status :: 3 - Alpha
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS
"""

params = {
    'name':             init.__name__,
    'version':          init.__version__,
    'license':          init.__license__,
    'long_description': open('README.md').read(),
    'long_description_content_type': 'text/markdown',
    'author':           init.__author__,
    'author_email':     init.__email__,
    'url':              init.__url__,
    'install_requires': ['biom-format'],
    'classifiers':      classes.strip().split('\n    '),
    'python_requires':  '>=3.6',
    'entry_points': {
        'console_scripts': [f'{init.__name__}=woltka.cli:cli'],
        'qiime2.plugins': [f'q2-{init.__name__}=woltka.q2.plugin_setup:plugin']
    }
}

setup(**params, packages=find_packages(), include_package_data=True)
