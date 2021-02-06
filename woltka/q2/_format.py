#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import model


class SeqAlnMapFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


SeqAlnMapDirFmt = model.SingleFileDirectoryFormat(
    'SeqAlnMapDirFmt', 'alignment.sam', SeqAlnMapFormat)


class BLAST6OutFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


BLAST6OutDirFmt = model.SingleFileDirectoryFormat(
    'BLAST6OutDirFmt', 'alignment.b6o', BLAST6OutFormat)


class SimpleMapFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


SimpleMapDirFmt = model.SingleFileDirectoryFormat(
    'SimpleMapDirFmt', 'mapping.txt', SimpleMapFormat)


class NCBINodesFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


NCBINodesDirFmt = model.SingleFileDirectoryFormat(
    'NCBINodesDirFmt', 'nodes.dmp', NCBINodesFormat)


class GeneCoordFormat(model.TextFileFormat):
    def _validate_(self, level):
        pass


GeneCoordDirFmt = model.SingleFileDirectoryFormat(
    'GeneCoordDirFmt', 'coordinates.txt', GeneCoordFormat)


# class NCBINamesFormat(model.TextFileFormat):
#     def _validate_(self, level):
#         pass


# class NCBITaxdumpDirFmt(model.DirectoryFormat):
#     nodes = model.File('nodes.dmp', format=NCBINodesFormat)
#     names = model.File('names.dmp', format=NCBINamesFormat)
