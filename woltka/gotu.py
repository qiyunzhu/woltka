#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import click

from woltka.classify import classify


@click.command()
@click.option(
    '--input', '-i', 'input_fp', required=True,
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    help='input read alignment directory')
@click.option(
    '--output', '-o', 'output_fp', required=True,
    type=click.Path(writable=True),
    help=('output gOTU table)'))
@click.option(
    '--multi/--no-multi', default=True,
    help=('allow one sequence to be assigned to multiple gOTUs; each hit '
          'will be counted as 1 / k (k is the totally number of hits)'))
@click.option(
    '--ixend/--no-ixend', default=False,
    help=('subject identifiers end with underscore index, the latter of which '
          'is to be removed prior to mapping.'))
@click.pass_context
def gotu(ctx, **kwargs):
    """Generate a gOTU table based on sequence alignments.
    """
    ctx.invoke(classify, **kwargs)


if __name__ == "__main__":
    gotu()
