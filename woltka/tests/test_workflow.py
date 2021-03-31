#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os import remove
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp
from filecmp import cmp

import pandas as pd
from biom import load_table
from pandas.testing import assert_frame_equal

from woltka.workflow import (
    workflow, classify, parse_samples, parse_strata, build_mapper,
    prepare_ranks, build_hierarchy, assign_readmap, strip_suffix, demultiplex,
    read_strata, round_profiles, write_profiles)


class WorkflowTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_workflow(self):
        # simplest ogu workflow
        input_fp = join(self.datdir, 'align', 'bowtie2')
        output_fp = join(self.tmpdir, 'tmp.tsv')
        obs = workflow(input_fp, output_fp)['none']
        self.assertEqual(obs['S01']['G000011545'], 48)
        self.assertNotIn('G000007145', obs['S02'])
        self.assertEqual(obs['S03']['G000009345'], 640)
        self.assertTrue(cmp(output_fp, join(
            self.datdir, 'output', 'bowtie2.ogu.tsv')))
        remove(output_fp)

    def test_classify(self):
        # simplest ogu workflow
        input_fp = join(self.datdir, 'align', 'bowtie2')
        samples, files, demux = parse_samples(input_fp)
        mapper, chunk = build_mapper()
        ranks = ['none']
        obs = classify(mapper, files, samples=samples, demux=demux,
                       ranks=ranks, chunk=chunk)['none']
        self.assertEqual(obs['S01']['G000011545'], 48)
        self.assertNotIn('G000007145', obs['S02'])
        self.assertEqual(obs['S03']['G000009345'], 640)
        self.assertEqual(obs['S04']['G000240185'], 4)
        self.assertEqual(obs['S05']['G000191145'], 10)

        # complex genus/process stratification workflow
        input_fp = join(self.datdir, 'align', 'burst')
        coords_fp = join(self.datdir, 'function', 'coords.txt.xz')
        map_fps = [join(self.datdir, 'function', 'uniref', 'uniref.map.xz'),
                   join(self.datdir, 'function', 'go', 'process.tsv.xz')]
        strata_dir = join(self.datdir, 'output', 'burst.genus.map')
        samples, files, demux = parse_samples(input_fp)
        tree, rankdic, namedic, root = build_hierarchy(
            map_fps=map_fps, map_as_rank=True)
        mapper, chunk = build_mapper(coords_fp=coords_fp, overlap=80)
        stratmap = parse_strata(strata_dir, samples)
        obs = classify(
            mapper, files, samples=samples, demux=demux, tree=tree,
            rankdic=rankdic, namedic=namedic, root=root, stratmap=stratmap,
            chunk=chunk, ranks=['process'])['process']
        self.assertEqual(obs['S01'][('Thermus', 'GO:0005978')], 2)
        self.assertEqual(obs['S02'][('Bacteroides', 'GO:0006814')], 1)
        self.assertEqual(obs['S03'][('Escherichia', 'GO:0006813')], 2)
        self.assertEqual(len(obs['S04']), 39)

    def test_parse_samples(self):
        # file (assuming demultiplexed)
        fp = join(self.tmpdir, 'input.fq')
        open(fp, 'a').close()
        obs = parse_samples(fp)
        self.assertIsNone(obs[0])
        self.assertListEqual(obs[1], [fp])
        self.assertTrue(obs[2])

        # confirm demux
        obs = parse_samples(fp, demux=True)
        self.assertIsNone(obs[0])
        self.assertListEqual(obs[1], [fp])
        self.assertTrue(obs[2])

        # force non-demux
        obs = parse_samples(fp, demux=False)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})
        self.assertFalse(obs[2])

        # provide correct extension
        obs = parse_samples(fp, ext='.fq', demux=False)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})

        # provide wrong extension
        with self.assertRaises(ValueError) as ctx:
            parse_samples(fp, ext='.fastq', demux=False)
        self.assertEqual(str(ctx.exception), (
            'Filepath and filename extension do not match.'))

        # provide Id list (no effect since demux)
        obs = parse_samples(fp, samples='input')
        self.assertListEqual(obs[0], ['input'])
        self.assertListEqual(obs[1], [fp])

        # provide correct Id list
        obs = parse_samples(fp, samples='input', demux=False)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})

        # provide wrong Id list
        with self.assertRaises(ValueError) as ctx:
            parse_samples(fp, samples='hello', demux=False)
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))

        # directory (assuming per-sample)
        obs = parse_samples(self.tmpdir)
        self.assertListEqual(obs[0], ['input'])
        self.assertDictEqual(obs[1], {fp: 'input'})
        self.assertFalse(obs[2])

        # directory empty
        remove(fp)
        with self.assertRaises(ValueError) as ctx:
            parse_samples(self.tmpdir)
        self.assertEqual(str(ctx.exception), (
            'No valid file found in directory.'))

        # multiple files
        for i in range(1, 4):
            open(join(self.tmpdir, f'S{i}.sam'), 'a').close()
        obs = parse_samples(self.tmpdir)[1]
        exp = {join(self.tmpdir, f'S{i}.sam'): f'S{i}' for i in range(1, 4)}
        self.assertDictEqual(obs, exp)

        # add an irrelevant file
        fp = join(self.tmpdir, 'readme.txt')
        open(fp, 'a').close()
        obs = parse_samples(self.tmpdir)[1]
        exp[fp] = 'readme'
        self.assertDictEqual(obs, exp)

        # specify extension to target alignment files only
        del exp[fp]
        obs = parse_samples(self.tmpdir, ext='.sam')[1]
        self.assertDictEqual(obs, exp)

        # specify sample Ids
        obs = parse_samples(self.tmpdir, samples='S1,S2,S3')[1]
        self.assertDictEqual(obs, exp)

        # some samples are not found
        remove(fp)
        with self.assertRaises(ValueError) as ctx:
            parse_samples(self.tmpdir, samples='S1,S2,S4')
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))

        # force demux
        obs = parse_samples(self.tmpdir, demux=True)
        self.assertIsNone(obs[0])
        self.assertTrue(obs[2])
        exp = [join(self.tmpdir, f'S{i}.sam') for i in range(1, 4)]
        self.assertListEqual(obs[1], exp)

        # sample Ids are ignored when demux
        obs = parse_samples(self.tmpdir, samples='S1,S2,S4', demux=True)
        self.assertListEqual(obs[0], ['S1', 'S2', 'S4'])
        exp = [join(self.tmpdir, f'S{i}.sam') for i in range(1, 4)]
        self.assertListEqual(obs[1], exp)

        # sample-to-file map
        fp = join(self.tmpdir, 'sample.list')
        with open(fp, 'w') as f:
            for i in (range(1, 4)):
                print(f'S{i}\tS{i}.sam', file=f)
        obs = parse_samples(fp)
        self.assertListEqual(obs[0], ['S1', 'S2', 'S3'])
        exp = {join(self.tmpdir, f'S{i}.sam'): f'S{i}' for i in range(1, 4)}
        self.assertDictEqual(obs[1], exp)

        # some samples only
        obs = parse_samples(fp, samples='S1,S2')
        self.assertListEqual(obs[0], ['S1', 'S2'])

        # some samples are not found
        with self.assertRaises(ValueError) as ctx:
            parse_samples(fp, samples='S1,S2,S4')
        self.assertEqual(str(ctx.exception), (
            'Provided sample IDs and actual files are inconsistent.'))
        remove(fp)

        # not a valid path
        with self.assertRaises(ValueError) as ctx:
            parse_samples('im/not/path')
        self.assertEqual(str(ctx.exception), (
            '"im/not/path" is not a valid file or directory.'))

    def test_parse_strata(self):
        # default
        for i in range(1, 4):
            open(join(self.tmpdir, f'S{i}.txt'), 'a').close()
        obs = parse_strata(self.tmpdir)
        exp = {f'S{i}': join(self.tmpdir, f'S{i}.txt') for i in range(1, 4)}
        self.assertDictEqual(obs, exp)

        # with sample Ids
        obs = parse_strata(self.tmpdir, samples=['S1', 'S2'])
        exp = {f'S{i}': join(self.tmpdir, f'S{i}.txt') for i in range(1, 3)}
        self.assertDictEqual(obs, exp)

        # sample missing
        with self.assertRaises(ValueError) as ctx:
            parse_strata(self.tmpdir, samples=['S1', 'S2', 'Sx'])
        self.assertEqual(str(ctx.exception), (
            'Cannot locate stratification files for one or more samples.'))

    def test_build_mapper(self):
        # plain
        obs = build_mapper()
        self.assertEqual(obs[0].__name__, 'plain_mapper')
        self.assertEqual(obs[1], 1000)

        # ordinal
        fp = join(self.tmpdir, 'coords.txt')
        with open(fp, 'w') as f:
            f.write('>G1\n1\t10\t20\n2\t35\t50\n')
        obs = build_mapper(fp)
        self.assertEqual(obs[0].func.__name__, 'ordinal_mapper')
        self.assertEqual(obs[1], 1000000)

        # ordinal with overlap threshold and chunk size
        obs = build_mapper(fp, overlap=75, chunk=50000)
        self.assertEqual(obs[0].func.__name__, 'ordinal_mapper')
        self.assertEqual(obs[0].keywords['th'], 0.75)
        self.assertEqual(obs[1], 50000)
        remove(fp)

    def test_prepare_ranks(self):
        # no rank
        obs = prepare_ranks()
        self.assertListEqual(obs[0], ['none'])

        # single rank
        obs = prepare_ranks('genus')
        self.assertListEqual(obs[0], ['genus'])

        # multiple ranks
        obs = prepare_ranks('phylum,genus,species')
        self.assertListEqual(obs[0], ['phylum', 'genus', 'species'])

        # missing rank
        with self.assertRaises(ValueError) as ctx:
            prepare_ranks('kingdom,phylum,class', rankdic={
                1: 'domain', 2: 'phylum', 3: 'phylum', 4: 'phylum'})
        self.assertEqual(str(ctx.exception), (
            'Ranks class, kingdom are not found in classification system.'))

        # free rank when there is tree
        obs = prepare_ranks(tree={'a': 1})
        self.assertListEqual(obs[0], ['free'])

        # with output read map directory
        obs = prepare_ranks(outmap_dir=self.tmpdir)
        self.assertTupleEqual(obs, (['none'], {'none': self.tmpdir}))

        # make multiple output read map directories
        obs = prepare_ranks('ko,mo', self.tmpdir)
        self.assertListEqual(obs[0], ['ko', 'mo'])
        self.assertDictEqual(obs[1], {'ko': join(self.tmpdir, 'ko'),
                                      'mo': join(self.tmpdir, 'mo')})
        rmtree(join(self.tmpdir, 'ko'))
        rmtree(join(self.tmpdir, 'mo'))

    def test_build_hierarchy(self):
        # names
        fp = join(self.tmpdir, 'names.dmp')
        with open(fp, 'w') as f:
            f.write('a\tHomo\nb\tPan\nc\tChimp\n')
        obs = build_hierarchy(names_fps=[fp])
        self.assertDictEqual(obs[2], {'a': 'Homo', 'b': 'Pan', 'c': 'Chimp'})
        remove(fp)

        # nodes
        fp = join(self.tmpdir, 'nodes.dmp')
        with open(fp, 'w') as f:
            f.write('a\td\nb\te\nc\td\nd\te\n')
        obs = build_hierarchy(nodes_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertEqual(obs[3], 'e')

        # nodes with rank
        fp = join(self.tmpdir, 'nodes.dmp')
        with open(fp, 'w') as f:
            f.write('a\td\tlv3\nb\te\tlv3\nc\td\tlv3\nd\te\tlv2\ne\te\tlv1\n')
        obs = build_hierarchy(nodes_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertDictEqual(obs[1], {'a': 'lv3', 'b': 'lv3', 'c': 'lv3',
                                      'd': 'lv2', 'e': 'lv1'})
        self.assertEqual(obs[3], 'e')
        remove(fp)

        # newick
        fp = join(self.tmpdir, 'tree.nwk')
        with open(fp, 'w') as f:
            f.write('((a,c)d,b)e;')
        obs = build_hierarchy(newick_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertEqual(obs[3], 'e')
        remove(fp)

        # lineage
        fp = join(self.tmpdir, 'lineage.txt')
        with open(fp, 'w') as f:
            f.write('a\te;d\nb\te\nc\te;d\n')
        obs = build_hierarchy(lineage_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'e;d', 'b': 'e', 'c': 'e;d', 'e;d': 'e', 'e': 'e'})
        self.assertEqual(obs[3], 'e')

        # lineage with rank code
        with open(fp, 'w') as f:
            f.write('a\tk__e;p__d\nb\tk__e\nc\tk__e;p__d\n')
        obs = build_hierarchy(lineage_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'k__e;p__d', 'b': 'k__e', 'c': 'k__e;p__d',
            'k__e;p__d': 'k__e', 'k__e': 'k__e'})
        self.assertDictEqual(obs[1], {
            'k__e': 'kingdom', 'k__e;p__d': 'phylum'})
        self.assertEqual(obs[3], 'k__e')
        remove(fp)

        # rank-per-column table
        fp = join(self.tmpdir, 'columns.tsv')
        with open(fp, 'w') as f:
            f.write('#ID\tlv1\tlv2\n'
                    'a\te\td\n'
                    'b\te\t\n'
                    'c\te\td\n')
        obs = build_hierarchy(columns_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertDictEqual(obs[1], {'d': 'lv2', 'e': 'lv1'})
        self.assertEqual(obs[3], 'e')
        remove(fp)

        # simple map
        fp = join(self.tmpdir, 'map.txt')
        with open(fp, 'w') as f:
            f.write('a\tBac\nb\tArc\nc\tBac\n')
        obs = build_hierarchy(map_fps=[fp])
        self.assertDictEqual(obs[0], {
            'a': 'Bac', 'b': 'Arc', 'c': 'Bac', 'Bac': '1', 'Arc': '1',
            '1': '1'})
        self.assertEqual(obs[3], '1')

        # map as rank
        obs = build_hierarchy(map_fps=[fp], map_as_rank=True)
        self.assertDictEqual(obs[1], {'Bac': 'map', 'Arc': 'map'})
        remove(fp)

    def test_strip_suffix(self):
        subs = [{'G1_1', 'G1_2', 'G2_3', 'G3'},
                {'G1_1', 'G1.3', 'G4_5', 'G4_x'}]
        obs = strip_suffix(subs, sep='_')
        exp = [{'G1', 'G2', 'G3'},
               {'G1', 'G1.3', 'G4', 'G4'}]
        self.assertListEqual(list(obs), exp)

        subs = [{'NC_123456.1_300', 'ABCD000001.20_101'}]
        obs = strip_suffix(subs, sep='_')
        exp = [{'NC_123456.1', 'ABCD000001.20'}]
        self.assertListEqual(list(obs), exp)

        subs = [{'G1.1', 'G1.2', 'G2'},
                {'G1.1', 'G1.3', 'G3_x'}]
        obs = strip_suffix(subs, sep='.')
        exp = [{'G1', 'G2'},
               {'G1', 'G3_x'}]
        self.assertListEqual(list(obs), exp)

    def test_demultiplex(self):
        # simple case
        rmap = [('S1_R1', 5),
                ('S1_R2', 12),
                ('S1_R3', 3),
                ('S2_R1', 10),
                ('S2_R2', 8),
                ('S2_R4', 7),
                ('S3_R2', 15),
                ('S3_R3', 1),
                ('S3_R4', 5)]
        obs = demultiplex(*zip(*rmap))
        exp = {'S1': [('R1',  5),
                      ('R2', 12),
                      ('R3',  3)],
               'S2': [('R1', 10),
                      ('R2',  8),
                      ('R4',  7)],
               'S3': [('R2', 15),
                      ('R3',  1),
                      ('R4',  5)]}
        self.assertEqual(obs.keys(), exp.keys())
        for s in obs:
            self.assertListEqual(list(map(tuple, obs[s])), list(zip(*exp[s])))

        # change separator, no result
        obs = demultiplex(*zip(*rmap), sep='.')
        self.assertEqual(obs.keys(), {''})
        self.assertListEqual(list(map(tuple, obs[''])), list(zip(*rmap)))

        # enforce sample Ids
        obs = demultiplex(*zip(*rmap), samples=['S1', 'S2', 'SX'])
        self.assertEqual(obs.keys(), {'S1', 'S2'})
        for s in ('S1', 'S2'):
            self.assertListEqual(list(map(tuple, obs[s])), list(zip(*exp[s])))

    def test_assign_readmap(self):
        # simple ogu assignment
        qryq = ['R1', 'R2', 'R3']
        subq = [frozenset(x) for x in [{'G1'}, {'G1', 'G2'}, {'G2', 'G3'}]]
        assigners = {}
        data = {'none': {}}
        assign_readmap(qryq, subq, data, 'none', 'S1', assigners)
        self.assertDictEqual(data['none']['S1'], {
            'G1': 1.5, 'G2': 1.0, 'G3': 0.5})

        # write read map
        data = {'none': {'S1': {}}}
        assign_readmap(qryq, subq, data, 'none', 'S1', assigners, rank2dir={
            'none': self.tmpdir})
        self.assertDictEqual(data['none']['S1'], {
            'G1': 1.5, 'G2': 1.0, 'G3': 0.5})
        fp = join(self.tmpdir, 'S1.txt')
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['R1\tG1', 'R2\tG1:1\tG2:1', 'R3\tG2:1\tG3:1']
        self.assertListEqual(obs, exp)
        remove(fp)

        # unique and unassigned
        assigners = {}
        data = {'none': {}}
        assign_readmap(qryq, subq, data, 'none', 'S1', assigners, uniq=True,
                       unasgd=True)
        self.assertDictEqual(data['none']['S1'], {'G1': 1, 'Unassigned': 2})

        # free-rank assignment
        tree = {'G1': 'T1', 'G2': 'T1', 'G3': 'T2',
                'T1': 'T0', 'T2': 'T0', 'T0': 'T0'}
        data = {'free': {}}
        assign_readmap(qryq, subq, data, 'free', 'S1', assigners, tree=tree)
        self.assertDictEqual(data['free']['S1'], {'T0': 1, 'T1': 2})

        # fixed-rank assignment
        rankdic = {'T1': 'ko', 'T2': 'ko', 'T0': 'mo'}
        data = {'ko': {}}
        assign_readmap(qryq, subq, data, 'ko', 'S1', assigners, tree=tree,
                       rankdic=rankdic)
        self.assertDictEqual(data['ko']['S1'], {'T1': 2.5, 'T2': 0.5})

    def test_read_strata(self):
        # regular strata file
        fp = join(self.datdir, 'output', 'burst.genus.map', 'S01.txt.gz')
        obs = read_strata(fp, zippers={'gzip': False})
        self.assertEqual(len(obs), 1608)

        # file has no strata
        with self.assertRaises(ValueError) as ctx:
            read_strata(join(self.datdir, 'tree.nwk'))
        self.assertEqual(str(ctx.exception), (
            'No stratification information is found in file: tree.nwk.'))

    def test_round_profiles(self):
        # free-rank: don't round
        obs = {'free': {'S1': {'G1': 1, 'G2': 2},
                        'S2': {'G1': 3, 'G3': 2}}}
        round_profiles(obs)
        exp = {'free': {'S1': {'G1': 1, 'G2': 2},
                        'S2': {'G1': 3, 'G3': 2}}}
        self.assertEqual(obs, exp)

        # none-rank: round
        obs = {'none': {'S1': {'G1': 1.3, 'G2': 2.2},
                        'S2': {'G1': 3.0, 'G3': 2.1}}}
        round_profiles(obs)
        exp = {'none': {'S1': {'G1': 1, 'G2': 2},
                        'S2': {'G1': 3, 'G3': 2}}}
        self.assertEqual(obs, exp)

        # none-rank and unique: don't round
        obs = {'none': {'S1': {'G1': 1.5, 'G2': 0},
                        'S2': {'G1': 3.0, 'G3': 2}}}
        round_profiles(obs, uniq=True)
        exp = {'none': {'S1': {'G1': 1.5, 'G2': 0},
                        'S2': {'G1': 3.0, 'G3': 2}}}
        self.assertDictEqual(obs, exp)

        # given-rank: round
        obs = {'class': {'S1': {'G1': 1.5, 'G2': 0},
                         'S2': {'G1': 3.0, 'G3': 2}}}
        round_profiles(obs)
        exp = {'class': {'S1': {'G1': 2},
                         'S2': {'G1': 3, 'G3': 2}}}
        self.assertDictEqual(obs, exp)

        # given-rank: don't round
        obs = {'class': {'S1': {'G1': 1.5, 'G2': 0},
                         'S2': {'G1': 3.0, 'G3': 2}}}
        round_profiles(obs, uniq=False, major=0, above=True)
        exp = {'class': {'S1': {'G1': 1.5, 'G2': 0},
                         'S2': {'G1': 3.0, 'G3': 2}}}
        self.assertDictEqual(obs, exp)

    def test_write_profiles(self):
        # do nothing
        self.assertIsNone(write_profiles({}, None))

        # single profile in TSV format
        fp = join(self.tmpdir, 'output.tsv')
        data = {'none': {'S1': {'G1': 1, 'G2': 2},
                         'S2': {'G1': 3, 'G3': 2}}}
        write_profiles(data, fp)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2', 'G1\t1\t3', 'G2\t2\t0', 'G3\t0\t2']
        self.assertListEqual(obs, exp)
        remove(fp)

        # force BIOM format (even if filename ends with .tsv)
        write_profiles(data, fp, is_biom=True)
        obs = load_table(fp).to_dataframe(dense=True).astype(int)
        exp = pd.DataFrame([[1, 3], [2, 0], [0, 2]], index=['G1', 'G2', 'G3'],
                           columns=['S1', 'S2'])
        assert_frame_equal(obs, exp)

        # multiple profiles in BIOM format (default)
        data = {'none': {'S1': {'G1': 1, 'G2': 2},
                         'S2': {'G1': 3, 'G3': 2}},
                'free': {'S1': {'Ecoli': 3, 'Strep': 0},
                         'S2': {'Ecoli': 1, 'Strep': 2}}}
        write_profiles(data, self.tmpdir)
        fp = join(self.tmpdir, 'none.biom')
        obs = load_table(fp).to_dataframe(dense=True).astype(int)
        assert_frame_equal(obs, exp)
        remove(fp)
        fp = join(self.tmpdir, 'free.biom')
        obs = load_table(fp).to_dataframe(dense=True).astype(int)
        exp = pd.DataFrame([[3, 1], [0, 2]], index=['Ecoli', 'Strep'],
                           columns=['S1', 'S2'])
        assert_frame_equal(obs, exp)
        remove(fp)


if __name__ == '__main__':
    main()
