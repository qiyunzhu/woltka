#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from .plugin_setup import plugin
from ._format import (SeqAlnMapDirFmt,
                      BLAST6OutDirFmt,
                      SimpleMapDirFmt,
                      NCBINodesDirFmt,
                      GeneCoordDirFmt)


@plugin.register_transformer
def _1(ff: SeqAlnMapDirFmt) -> str:
    return join(str(ff.path), 'alignment.sam')


@plugin.register_transformer
def _2(ff: BLAST6OutDirFmt) -> str:
    return join(str(ff.path), 'alignment.b6o')


@plugin.register_transformer
def _3(ff: SimpleMapDirFmt) -> str:
    return join(str(ff.path), 'mapping.txt')


@plugin.register_transformer
def _4(ff: NCBINodesDirFmt) -> str:
    return join(str(ff.path), 'nodes.dmp')


@plugin.register_transformer
def _5(ff: GeneCoordDirFmt) -> str:
    return join(str(ff.path), 'coordinates.txt')
