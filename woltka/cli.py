#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click
from importlib import import_module


click_kws = {
    'context_settings': {
        'help_option_names': ['-h', '--help']}}

modules = ('classify', 'gotu')


@click.version_option('0.1.0')
@click.group(**click_kws)
def main():
    pass


for module in modules:
    main.add_command(getattr(import_module(f'woltka.{module}'), module))


if __name__ == '__main__':
    main()
