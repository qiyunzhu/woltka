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

from woltka.workflow import (
    workflow, classify, parse_samples, parse_strata, build_mapper,
    prepare_ranks, build_hierarchy, reshape_readmap, assign_readmap,
    write_profiles)


class WorkflowTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_workflow(self):
        # simplest gotu workflow
        input_path = join(self.datdir, 'align', 'bowtie2')
        output_path = join(self.tmpdir, 'tmp.tsv')
        obs = workflow(input_path, output_path)['none']
        self.assertEqual(obs['S01']['G000011545'], 48)
        self.assertNotIn('G000007145', obs['S02'])
        self.assertEqual(obs['S03']['G000009345'], 640)
        with open(output_path, 'r') as f:
            obs = f.read().splitlines()
        exp_fp = join(self.datdir, 'output', 'bowtie2.gotu.tsv')
        with open(exp_fp, 'r') as f:
            exp = f.read().splitlines()
        self.assertListEqual(obs, exp)
        remove(output_path)

    def test_classify(self):
        # simplest gotu workflow
        input_path = join(self.datdir, 'align', 'bowtie2')
        samples, files, demux = parse_samples(input_path)
        mapper = build_mapper()
        ranks = ['none']
        obs = classify(mapper, files, samples=samples, demux=demux,
                       ranks=ranks)['none']
        self.assertEqual(obs['S01']['G000011545'], 48)
        self.assertNotIn('G000007145', obs['S02'])
        self.assertEqual(obs['S03']['G000009345'], 640)
        self.assertEqual(obs['S04']['G000240185'], 4)
        self.assertEqual(obs['S05']['G000191145'], 10)

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
        self.assertEqual(obs.__class__.__name__, 'Plain')

        # ordinal
        fp = join(self.tmpdir, 'coords.txt')
        with open(fp, 'w') as f:
            f.write('>G1\n1\t10\t20\n2\t35\t50\n')
        obs = build_mapper(fp)
        self.assertEqual(obs.__class__.__name__, 'Ordinal')

        # ordinal with overlap threshold
        obs = build_mapper(fp, 75)
        self.assertEqual(obs.__class__.__name__, 'Ordinal')
        self.assertEqual(obs.th, 0.75)
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
        obs = build_hierarchy(names_fp=fp)
        self.assertDictEqual(obs[2], {'a': 'Homo', 'b': 'Pan', 'c': 'Chimp'})
        remove(fp)

        # nodes
        fp = join(self.tmpdir, 'nodes.dmp')
        with open(fp, 'w') as f:
            f.write('a\td\nb\te\nc\td\nd\te\n')
        obs = build_hierarchy(nodes_fp=fp)
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertEqual(obs[3], 'e')

        # nodes with rank
        fp = join(self.tmpdir, 'nodes.dmp')
        with open(fp, 'w') as f:
            f.write('a\td\tlv3\nb\te\tlv3\nc\td\tlv3\nd\te\tlv2\ne\te\tlv1\n')
        obs = build_hierarchy(nodes_fp=fp)
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
        obs = build_hierarchy(newick_fp=fp)
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertEqual(obs[3], 'e')
        remove(fp)

        # lineage
        fp = join(self.tmpdir, 'lineage.txt')
        with open(fp, 'w') as f:
            f.write('a\te;d\nb\te\nc\te;d\n')
        obs = build_hierarchy(lineage_fp=fp)
        self.assertDictEqual(obs[0], {
            'a': 'e;d', 'b': 'e', 'c': 'e;d', 'e;d': 'e', 'e': 'e'})
        self.assertEqual(obs[3], 'e')

        # lineage with rank code
        with open(fp, 'w') as f:
            f.write('a\tk__e;p__d\nb\tk__e\nc\tk__e;p__d\n')
        obs = build_hierarchy(lineage_fp=fp)
        self.assertDictEqual(obs[0], {
            'a': 'k__e;p__d', 'b': 'k__e', 'c': 'k__e;p__d',
            'k__e;p__d': 'k__e', 'k__e': 'k__e'})
        self.assertDictEqual(obs[1], {
            'k__e': 'kingdom', 'k__e;p__d': 'phylum'})
        self.assertEqual(obs[3], 'k__e')
        remove(fp)

        # rank table
        fp = join(self.tmpdir, 'ranks.tsv')
        with open(fp, 'w') as f:
            f.write('#ID\tlv1\tlv2\n'
                    'a\te\td\n'
                    'b\te\t\n'
                    'c\te\td\n')
        obs = build_hierarchy(ranktb_fp=fp)
        self.assertDictEqual(obs[0], {
            'a': 'd', 'b': 'e', 'c': 'd', 'd': 'e', 'e': 'e'})
        self.assertDictEqual(obs[1], {'d': 'lv2', 'e': 'lv1'})
        self.assertEqual(obs[3], 'e')
        remove(fp)

        # map
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

    def test_reshape_readmap(self):
        # doing nothing
        rmap = {'R1': {'G1'}, 'R2': {'G2'}, 'R3': {'G3', 'G4'}}
        obs = reshape_readmap(rmap)
        self.assertDictEqual(obs, {None: rmap})

        # with filename
        obs = reshape_readmap(rmap, files={'fname': 'S1'}, fp='fname')
        self.assertDictEqual(obs, {'S1': rmap})

        # remove index
        rmap_ = {'R1': {'G1'}, 'R2': {'G2_1'}, 'R3': {'G3_1', 'G4_2'}}
        obs = reshape_readmap(rmap_, deidx=True)
        self.assertDictEqual(obs, {None: rmap})

        # demultiplex
        rmap = {'S1_R1': {'G1'}, 'S1_R2': {'G2'}, 'S2_R3': {'G3', 'G4'}}
        obs = reshape_readmap(rmap, demux=True)
        self.assertDictEqual(obs, {'S1': {'R1': {'G1'}, 'R2': {'G2'}},
                                   'S2': {'R3': {'G3', 'G4'}}})

        # demultiplex to given sample Ids
        rmap = {'S1_R1': {'G1'}, 'S1_R2': {'G2'}, 'S2_R3': {'G3', 'G4'}}
        obs = reshape_readmap(rmap, demux=True, samples=['S1'])
        self.assertDictEqual(obs, {'S1': {'R1': {'G1'}, 'R2': {'G2'}}})

    def test_assign_readmap(self):
        # simple gotu assignment
        rmap = {'R1': {'G1'}, 'R2': {'G1', 'G2'}, 'R3': {'G2', 'G3'}}
        data = {'none': {}}
        assign_readmap(rmap, data, 'none', 'S1')
        self.assertDictEqual(data['none']['S1'], {'G1': 2, 'G2': 1})

        # write read map
        data = {'none': {'S1': {}}}
        assign_readmap(rmap, data, 'none', 'S1', rank2dir={
            'none': self.tmpdir})
        self.assertDictEqual(data['none']['S1'], {'G1': 2, 'G2': 1})
        fp = join(self.tmpdir, 'S1.txt')
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['R1\tG1', 'R2\tG1:1\tG2:1', 'R3\tG2:1\tG3:1']
        self.assertListEqual(obs, exp)
        remove(fp)

        # free-rank assignment
        tree = {'G1': 'T1', 'G2': 'T1', 'G3': 'T2',
                'T1': 'T0', 'T2': 'T0', 'T0': 'T0'}
        data = {'free': {}}
        assign_readmap(rmap, data, 'free', 'S1', tree=tree)
        self.assertDictEqual(data['free']['S1'], {'T0': 1, 'T1': 2})

        # fixed-rank assignment
        rankdic = {'T1': 'ko', 'T2': 'ko', 'T0': 'mo'}
        data = {'ko': {}}
        assign_readmap(rmap, data, 'ko', 'S1', tree=tree, rankdic=rankdic)
        self.assertDictEqual(data['ko']['S1'], {'T1': 2})

    def test_write_profiles(self):
        # do nothing
        self.assertIsNone(write_profiles({}, None))

        # simple
        fp = join(self.tmpdir, 'output.tsv')
        data = {'none': {'S1': {'G1': 1, 'G2': 2},
                         'S2': {'G1': 3, 'G3': 2}}}
        write_profiles(data, fp)
        with open(fp, 'r') as f:
            obs = f.read().splitlines()
        exp = ['#FeatureID\tS1\tS2', 'G1\t1\t3', 'G2\t2\t0', 'G3\t0\t2']
        self.assertListEqual(obs, exp)
        remove(fp)


if __name__ == '__main__':
    main()
