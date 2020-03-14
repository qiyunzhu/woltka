#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import setup, find_packages
from re import search


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

classifiers = [s.strip() for s in classes.split('\n') if s]


name = 'woltka'

description = 'Web of Life ToolKit App'

long_description = open('README.md').read()

with open(f'{name}/__init__.py', 'r') as f:
    version = search(r'__version__ = (.*)', f.read()).group(1).strip('\'"')

setup(
    name=name,
    version=version,
    license='BSD-3-Clause',
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Qiyun Zhu',
    author_email='qiyunzhu@gmail.com',
    url=f'https://github.com/qiyunzhu/{name}',
    packages=find_packages(),
    package_data={name: ['config.yml']},
    include_package_data=True,
    install_requires=[
        'biom-format',
        'cython'
    ],
    classifiers=classifiers,
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [f'{name}={name}.cli:cli'],
        'qiime2.plugins': [f'q2-{name}={name}.q2.plugin_setup:plugin']
    }
)
