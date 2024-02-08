#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin import model, ValidationError

from woltka.align import (
    check_map_file, check_b6o_file, check_sam_file, check_paf_file)
from woltka.ordinal import check_gene_coords


class SeqAlnMapFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.path.open('r') as fh:
            try:
                check_sam_file(fh)
            except ValueError as err:
                raise ValidationError(str(err))


SeqAlnMapDirFmt = model.SingleFileDirectoryFormat(
    'SeqAlnMapDirFmt', 'alignment.sam', SeqAlnMapFormat)


class PaimApFmtFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.path.open('r') as fh:
            try:
                check_paf_file(fh)
            except ValueError as err:
                raise ValidationError(str(err))


PaimApFmtDirFmt = model.SingleFileDirectoryFormat(
    'PaimApFmtDirFmt', 'alignment.paf', PaimApFmtFormat)


class BLAST6OutFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.path.open('r') as fh:
            try:
                check_b6o_file(fh)
            except ValueError as err:
                raise ValidationError(str(err))


BLAST6OutDirFmt = model.SingleFileDirectoryFormat(
    'BLAST6OutDirFmt', 'alignment.b6o', BLAST6OutFormat)


class SimpleMapFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.path.open('r') as fh:
            try:
                check_map_file(fh)
            except ValueError as err:
                raise ValidationError(str(err))


SimpleMapDirFmt = model.SingleFileDirectoryFormat(
    'SimpleMapDirFmt', 'mapping.txt', SimpleMapFormat)


class NCBINodesFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.path.open('r') as fh:
            for line in fh:
                line = line.rstrip()
                x = line.replace('\t|', '').split('\t')
                if len(x) < 2 or not x[0] or not x[1]:
                    raise ValidationError('Invalid node mapping: "{line}".')


NCBINodesDirFmt = model.SingleFileDirectoryFormat(
    'NCBINodesDirFmt', 'nodes.dmp', NCBINodesFormat)


class GeneCoordFormat(model.TextFileFormat):
    def _validate_(self, level):
        with self.path.open('r') as fh:
            try:
                check_gene_coords(fh)
            except ValueError as err:
                raise ValidationError(str(err))


GeneCoordDirFmt = model.SingleFileDirectoryFormat(
    'GeneCoordDirFmt', 'coordinates.txt', GeneCoordFormat)


# class NCBINamesFormat(model.TextFileFormat):
#     def _validate_(self, level):
#         pass


# class NCBITaxdumpDirFmt(model.DirectoryFormat):
#     nodes = model.File('nodes.dmp', format=NCBINodesFormat)
#     names = model.File('names.dmp', format=NCBINamesFormat)
