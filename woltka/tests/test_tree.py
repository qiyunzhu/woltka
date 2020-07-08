#!/usr/bin/env python3

# ----------------------------------------------------------------------------
# Copyright (c) 2020--, Qiyun Zhu.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from unittest import TestCase, main
from os.path import join, dirname, realpath
from shutil import rmtree
from tempfile import mkdtemp

from woltka.tree import (
    read_names, read_nodes, read_lineage, read_newick, read_columns,
    fill_root, get_lineage, lineage_str, find_rank, find_lca)


# A small test taxon set containing 15 proteobacterial species

PROTEO = {'taxmap': {  # genome Id: taxonomy Id
    'G000008865': '562',    # Escherichia coli
    'G000195995': '28901',  # Salmonella enterica
    'G000240185': '573',    # Klebsiella pneumoniae
    'G000025565': '550',    # Enterobacter cloacae
    'G000009065': '632',    # Yersinia pestis
    'G000006745': '666',    # Vibrio cholerae
    'G000012245': '317',    # Pseudomonas syringae
    'G000191145': '48296',  # Acinetobacter pittii
    'G000195715': '520',    # Bordetella pertussis
    'G000011705': '13373',  # Burkholderia mallei
    'G000006965': '382',    # Sinorhizobium meliloti
    'G000012905': '1063',   # Rhodobacter sphaeroides
    'G000195735': '782',    # Rickettsia prowazekii
    'G000008525': '210',    # Helicobacter pylori
    'G000195755': '881'},   # Desulfovibrio vulgaris

          'phylogeny': (  # Newick format
    '((((((G000191145:1.436,G000012245:0.929)N63:0.189,((((G000195995:0.140,G0'
    '00008865:0.097)N94:0.033,(G000025565:0.150,G000240185:0.127)N95:0.027)N91'
    ':0.220,G000009065:0.266)N84:0.607,G000006745:0.737)N73:1.229)N53:0.814,(G'
    '000195715:0.707,G000011705:0.657)N46:0.835)N28:1.647,(G000195735:2.210,(G'
    '000006965:1.256,G000012905:1.058)N48:0.544)N29:1.044)N20:1.107,G000195755'
    ':2.108)N13:1.043,G000008525:3.548)N2:0.207;'),

          'taxonomy': (  # Id, parent, rank, name
    '209',     '72293',   'genus',   'Helicobacter',
    '210',     '209',     'species', 'Helicobacter pylori',
    '286',     '135621',  'genus',   'Pseudomonas',
    '317',     '286',     'species', 'Pseudomonas syringae',
    '356',     '28211',   'order',   'Rhizobiales',
    '382',     '28105',   'species', 'Sinorhizobium meliloti',
    '468',     '72274',   'family',  'Moraxellaceae',
    '469',     '468',     'genus',   'Acinetobacter',
    '506',     '80840',   'family',  'Alcaligenaceae',
    '517',     '506',     'genus',   'Bordetella',
    '520',     '517',     'species', 'Bordetella pertussis',
    '543',     '91347',   'family',  'Enterobacteriaceae',
    '547',     '543',     'genus',   'Enterobacter',
    '550',     '547',     'species', 'Enterobacter cloacae',
    '561',     '543',     'genus',   'Escherichia',
    '562',     '561',     'species', 'Escherichia coli',
    '570',     '543',     'genus',   'Klebsiella',
    '573',     '570',     'species', 'Klebsiella pneumoniae',
    '590',     '543',     'genus',   'Salmonella',
    '629',     '1903411', 'genus',   'Yersinia',
    '632',     '629',     'species', 'Yersinia pestis',
    '641',     '135623',  'family',  'Vibrionaceae',
    '662',     '641',     'genus',   'Vibrio',
    '666',     '662',     'species', 'Vibrio cholerae',
    '766',     '28211',   'order',   'Rickettsiales',
    '775',     '766',     'family',  'Rickettsiaceae',
    '780',     '775',     'genus',   'Rickettsia',
    '782',     '780',     'species', 'Rickettsia prowazekii',
    '872',     '194924',  'genus',   'Desulfovibrio',
    '881',     '872',     'species', 'Desulfovibrio vulgaris',
    '1060',    '31989',   'genus',   'Rhodobacter',
    '1063',    '1060',    'species', 'Rhodobacter sphaeroides',
    '1224',    '1224',    'phylum',  'Proteobacteria',
    '1236',    '1224',    'class',   'Gammaproteobacteria',
    '13373',   '32008',   'species', 'Burkholderia mallei',
    '28105',   '82115',   'genus',   'Sinorhizobium',
    '28211',   '1224',    'class',   'Alphaproteobacteria',
    '28216',   '1224',    'class',   'Betaproteobacteria',
    '28221',   '1224',    'class',   'Deltaproteobacteria',
    '28901',   '590',     'species', 'Salmonella enterica',
    '29547',   '1224',    'class',   'Epsilonproteobacteria',
    '31989',   '204455',  'family',  'Rhodobacteraceae',
    '32008',   '119060',  'genus',   'Burkholderia',
    '48296',   '469',     'species', 'Acinetobacter pittii',
    '72274',   '1236',    'order',   'Pseudomonadales',
    '72293',   '213849',  'family',  'Helicobacteraceae',
    '80840',   '28216',   'order',   'Burkholderiales',
    '82115',   '356',     'family',  'Rhizobiaceae',
    '91347',   '1236',    'order',   'Enterobacterales',
    '119060',  '80840',   'family',  'Burkholderiaceae',
    '135621',  '72274',   'family',  'Pseudomonadaceae',
    '135623',  '1236',    'order',   'Vibrionales',
    '194924',  '213115',  'family',  'Desulfovibrionaceae',
    '204455',  '28211',   'order',   'Rhodobacterales',
    '213115',  '28221',   'order',   'Desulfovibrionales',
    '213849',  '29547',   'order',   'Campylobacterales',
    '1903411', '91347',   'family',  'Yersiniaceae')}

# tree diagram
#
#                                 /-G000191145 Acinetobacter
#                         /N63---|
#                        |        \-G000012245 Pseudomonas
#                        |
#                        |              /-G000195995 Salmonella
#                        |           /N94
#                     /N53          |   \-G000008865 Escherichia
#                    |   |        /N91
#                    |   |       |  |   /-G000025565 Enterobacter
#                    |   |    /N84   \N95
#                    |   |   |   |      \-G000240185 Klebsiella
#                 /N28    \N73   |
#                |   |       |    \-G000009065 Yersinia
#                |   |       |
#                |   |        \-----G000006745 Vibrio
#                |   |
#             /N20   |          /-G000195715 Bordetella
#            |   |    \N46-----|
#            |   |              \-G000011705 Burkholderia
#            |   |
#            |   |          /-G000195735 Rickettsia
#         /N13    \N29-----|
#        |   |             |        /-G000006965 Sinorhizobium
#        |   |              \N48---|
# -N2----|   |                      \-G000012905 Rhodobacter
#        |   |
#        |    \---------------G000195755 Desulfovibrio
#        |
#         \---------------G000008525 Helicobacter


class TreeTests(TestCase):
    def setUp(self):
        self.tmpdir = mkdtemp()
        self.datdir = join(dirname(realpath(__file__)), 'data')
        lst_ = PROTEO['taxonomy']
        ran_ = range(0, len(lst_), 4)
        self.proteo = {
            'map': PROTEO['taxmap'],
            'tree':   dict((lst_[i], lst_[i + 1]) for i in ran_),
            'ranks':  dict((lst_[i], lst_[i + 2]) for i in ran_),
            'names':  dict((lst_[i], lst_[i + 3]) for i in ran_)}

    def tearDown(self):
        rmtree(self.tmpdir)

    def test_read_names(self):
        # simple map
        tsv = ('1	root',
               '2	Bacteria',
               '6	Azorhizobium	* notes')
        obs = read_names(tsv)
        exp = {'1': 'root', '2': 'Bacteria', '6': 'Azorhizobium'}
        self.assertDictEqual(obs, exp)

        # names.dmp
        dmp = ('2	|	Bacteria	|	Bacteria	|	scientific name	|',
               '2	|	Procaryotae	|	Procaryotae <Bacteria>	|	in-part	|',
               '2157	|	Woese et al. 1990	|		|	authority	|',
               '2157	|	Archaea	|		|	scientific name	|')
        obs = read_names(dmp)
        exp = {'2': 'Bacteria', '2157': 'Archaea'}
        self.assertDictEqual(obs, exp)

        # real names file
        fp = join(self.datdir, 'taxonomy', 'names.dmp')
        with open(fp, 'r') as f:
            obs = read_names(f)
        self.assertEqual(len(obs), 486)
        self.assertEqual(obs['356'], 'Rhizobiales')
        self.assertEqual(obs['1117'], 'Cyanobacteria')

    def test_read_nodes(self):
        # simple map
        tsv = ('1	1',
               '2	1',
               '1224	2',
               '1236	1224',
               '32066	2')
        obs0 = read_nodes(tsv)[0]
        exp0 = {'1': '1', '2': '1', '1224': '2', '1236': '1224', '32066': '2'}
        self.assertDictEqual(obs0, exp0)

        # simple map with rank
        tsv = ('1	1	no rank',
               '2	1	superkingdom',
               '1224	2	phylum',
               '1236	1224	class',
               '32066	2	phylum')
        obs0, obs1 = read_nodes(tsv)
        self.assertDictEqual(obs0, exp0)
        exp1 = {'1': 'no rank', '2': 'superkingdom', '1224': 'phylum',
                '1236': 'class', '32066': 'phylum'}
        self.assertDictEqual(obs1, exp1)

        # nodes.dmp
        dmp = ('1	|	1	|	no rank	|',
               '2	|	1	|	superkingdom	|',
               '1224	|	2	|	phylum	|',
               '1236	|	1224	|	class	|',
               '32066	|	2	|	phylum	|')
        obs0, obs1 = read_nodes(dmp)
        self.assertDictEqual(obs0, exp0)
        self.assertDictEqual(obs1, exp1)

        # real nodes file
        fp = join(self.datdir, 'taxonomy', 'nodes.dmp')
        with open(fp, 'r') as f:
            obs0, obs1 = read_nodes(f)
        self.assertEqual(len(obs0), 486)
        self.assertEqual(obs0['356'], '28211')
        self.assertEqual(obs0['1117'], '1798711')
        self.assertEqual(len(obs1), 486)
        self.assertEqual(obs1['561'], 'genus')
        self.assertEqual(obs1['562'], 'species')

    def test_read_newick(self):
        # simple tree
        nwk = '((a,b)c,(d,e)f)g;'
        exp = {'a': 'c', 'b': 'c', 'd': 'f', 'e': 'f',
               'c': 'g', 'f': 'g', 'g': 'g'}
        self.assertDictEqual(read_newick((nwk,)), exp)

        # tree with branch lengths
        nwk = '((a:0.25,b:0.12)c:2.0,(d,e:1.75)f:0.55)g;'
        self.assertDictEqual(read_newick((nwk,)), exp)

        # multi-line tree file
        nwk = (' ((a:0.25,b:0.12)c:2.0, '
               ' (d,e:1.75)f:0.55)g;    ')
        self.assertDictEqual(read_newick(nwk), exp)

        # tree without internal node Id
        nwk = '((a,b),(c,d));'
        with self.assertRaises(ValueError) as ctx:
            read_newick((nwk,))
        self.assertEqual(str(ctx.exception), 'Missing internal node ID.')

        # tree with duplicate node Ids
        nwk = '((a,b)x,(c,d)x)y;'
        with self.assertRaises(ValueError) as ctx:
            read_newick((nwk,))
        self.assertEqual(str(ctx.exception), 'Found non-unique node ID: "x".')

        # proteo tree
        obs = read_newick((PROTEO['phylogeny'],))
        exp = {'G000006965': 'N48', 'G000012905': 'N48', 'G000195735': 'N29',
               'N48': 'N29', 'G000195715': 'N46', 'G000011705': 'N46',
               'G000025565': 'N95', 'G000240185': 'N95', 'G000195995': 'N94',
               'G000008865': 'N94', 'N94': 'N91', 'N95': 'N91', 'N91': 'N84',
               'G000009065': 'N84', 'N84': 'N73', 'G000006745': 'N73',
               'G000191145': 'N63', 'G000012245': 'N63', 'N63': 'N53', 'N73':
               'N53', 'N53': 'N28', 'N46': 'N28', 'N28': 'N20', 'N29': 'N20',
               'N20': 'N13', 'G000195755': 'N13', 'N13': 'N2', 'G000008525':
               'N2', 'N2': 'N2'}
        self.assertDictEqual(obs, exp)

        # real tree file
        fp = join(self.datdir, 'tree.nwk')
        with open(fp, 'r') as f:
            obs = read_newick(f)
        self.assertEqual(len(obs), 213)
        self.assertEqual(obs['G000027165'], 'N101')
        self.assertEqual(obs['N55'], 'N46')
        self.assertEqual(obs['N2'], 'N1')
        self.assertEqual(obs['N1'], 'N1')

    def test_read_columns(self):
        # simple case
        tsv = ['#seq	k	p	c	o',
               'seq1	k1	p1	c1	o1',
               'seq2	k1	p1	c2	o2',
               'seq3	k1	p2	c3	',
               'seq4	k2		c4	o3']
        tree_obs, rankdic_obs = read_columns(iter(tsv))
        tree_exp = {'seq1': 'o1', 'o1': 'c1', 'c1': 'p1', 'p1': 'k1',
                    'seq2': 'o2', 'o2': 'c2', 'c2': 'p1',
                    'seq3': 'c3', 'c3': 'p2', 'p2': 'k1',
                    'seq4': 'o3', 'o3': 'c4', 'c4': 'k2',
                    'k1': None, 'k2': None}
        rankdic_exp = {'k1': 'k', 'k2': 'k',
                       'p1': 'p', 'p2': 'p',
                       'c1': 'c', 'c2': 'c', 'c3': 'c', 'c4': 'c',
                       'o1': 'o', 'o2': 'o', 'o3': 'o'}
        self.assertDictEqual(tree_obs, tree_exp)
        self.assertDictEqual(rankdic_obs, rankdic_exp)

        # inconsistent tree
        with self.assertRaises(ValueError) as ctx:
            read_columns(iter(tsv + ['seq5	k1	p2	c1	o1']))
        self.assertEqual(str(ctx.exception), 'Conflict at taxon "c1".')

        # inconsistent rank
        with self.assertRaises(ValueError) as ctx:
            read_columns(iter(tsv + ['seq5	k2	c4		o3']))
        self.assertEqual(str(ctx.exception), 'Conflict at taxon "c4".')

        # real rank table file
        fp = join(self.datdir, 'taxonomy', 'rank_tids.tsv')
        with open(fp, 'r') as f:
            tree, rankdic = read_columns(f)
        self.assertEqual(len(tree), 426)
        self.assertEqual(len(rankdic), 319)
        self.assertEqual(tree['562'], '561')
        self.assertEqual(tree['356'], '28211')
        self.assertEqual(rankdic['543'], 'family')
        self.assertEqual(rankdic['1224'], 'phylum')
        self.assertIsNone(tree['2'])

    def test_read_lineage(self):
        # simple case
        tsv = ('seq1	k1;p1;c1;o1',
               'seq2	k1;p2;c2;o2',
               'seq3	k1;p2;c3;',    # empty taxon in end
               'seq4	k2;p3;c1;o3',  # duplicate taxon "c1"
               'seq5	k2;;c4;o4')    # empty taxon in middle
        obs = read_lineage(tsv)[0]
        exp = {'k1':          None,
               'k2':          None,
               'k1;p1':       'k1',
               'k1;p2':       'k1',
               'k2;p3':       'k2',
               'k1;p1;c1':    'k1;p1',
               'k1;p2;c2':    'k1;p2',
               'k1;p2;c3':    'k1;p2',
               'k2;p3;c1':    'k2;p3',
               'k2;;c4':      'k2',
               'k1;p1;c1;o1': 'k1;p1;c1',
               'k1;p2;c2;o2': 'k1;p2;c2',
               'k2;p3;c1;o3': 'k2;p3;c1',
               'k2;;c4;o4':   'k2;;c4',
               'seq1':        'k1;p1;c1;o1',
               'seq2':        'k1;p2;c2;o2',
               'seq3':        'k1;p2;c3',
               'seq4':        'k2;p3;c1;o3',
               'seq5':        'k2;;c4;o4'}
        self.assertDictEqual(obs, exp)

        # simple case with rank code
        tsv = ('#OTU ID	Taxon',
               'G1	k__Bacteria; p__Firmicutes; c__Bacilli',
               'G2	k__Bacteria; p__Firmicutes; c__Clostridia',
               'G3	k__Bacteria; p__Tenericutes; c__Mollicutes',
               'G4	k__Viruses; x__Inoviridae')
        tree_obs, rankdic_obs = read_lineage(tsv)
        tree_exp = {'k__Bacteria': None,
                    'k__Viruses': None,
                    'k__Bacteria;p__Firmicutes': 'k__Bacteria',
                    'k__Bacteria;p__Tenericutes': 'k__Bacteria',
                    'k__Bacteria;p__Firmicutes;c__Bacilli':
                    'k__Bacteria;p__Firmicutes',
                    'k__Bacteria;p__Firmicutes;c__Clostridia':
                    'k__Bacteria;p__Firmicutes',
                    'k__Bacteria;p__Tenericutes;c__Mollicutes':
                    'k__Bacteria;p__Tenericutes',
                    'k__Viruses;x__Inoviridae': 'k__Viruses',
                    'G1': 'k__Bacteria;p__Firmicutes;c__Bacilli',
                    'G2': 'k__Bacteria;p__Firmicutes;c__Clostridia',
                    'G3': 'k__Bacteria;p__Tenericutes;c__Mollicutes',
                    'G4': 'k__Viruses;x__Inoviridae'}
        self.assertDictEqual(tree_obs, tree_exp)
        rankdic_exp = {'k__Bacteria': 'kingdom',
                       'k__Viruses': 'kingdom',
                       'k__Bacteria;p__Firmicutes': 'phylum',
                       'k__Bacteria;p__Tenericutes': 'phylum',
                       'k__Bacteria;p__Firmicutes;c__Bacilli': 'class',
                       'k__Bacteria;p__Firmicutes;c__Clostridia': 'class',
                       'k__Bacteria;p__Tenericutes;c__Mollicutes': 'class'}
        self.assertDictEqual(rankdic_obs, rankdic_exp)

        # real lineage file
        fp = join(self.datdir, 'taxonomy', 'lineage.txt')
        with open(fp, 'r') as f:
            tree, rankdic = read_lineage(f)
        self.assertEqual(len(tree), 426)
        self.assertEqual(len(rankdic), 319)
        self.assertEqual(tree['k__Bacteria;p__Firmicutes'], 'k__Bacteria')
        self.assertEqual(rankdic['k__Bacteria;p__Firmicutes'], 'phylum')
        self.assertIsNone(tree['k__Bacteria'])

    def test_get_lineage(self):
        # simple tree
        tsv = ('1	1',
               '2	1',
               '1224	2',
               '1236	1224',
               '28211	1224'
               '32066	2')
        tree = read_nodes(tsv)[0]
        obs = get_lineage('1236', tree)
        exp = ['1', '2', '1224', '1236']
        self.assertListEqual(obs, exp)
        self.assertIsNone(get_lineage('1234', tree))

        # proteo tree
        obs = get_lineage('561', self.proteo['tree'])  # Escherichia
        exp = ['1224', '1236', '91347', '543', '561']
        self.assertListEqual(obs, exp)
        obs = get_lineage('286', self.proteo['tree'])  # Pseudomonas
        exp = ['1224', '1236', '72274', '135621', '286']
        self.assertListEqual(obs, exp)

        # real tree file
        fp = join(self.datdir, 'taxonomy', 'nodes.dmp')
        with open(fp, 'r') as f:
            tree = read_nodes(f)[0]
        obs = get_lineage('561', tree)  # Escherichia
        exp = ['1', '131567', '2', '1224', '1236', '91347', '543', '561']
        self.assertListEqual(obs, exp)
        obs = get_lineage('286', tree)  # Pseudomonas
        exp = ['1', '131567', '2', '1224', '1236', '72274', '135621', '286']
        self.assertListEqual(obs, exp)

    def test_lineage_str(self):
        tree = self.proteo['tree']
        obs = lineage_str('561', tree)
        self.assertEqual(obs, '1236;91347;543')
        obs = lineage_str('561', tree, include_root=True)
        self.assertEqual(obs, '1224;1236;91347;543')
        obs = lineage_str('561', tree, include_self=True)
        self.assertEqual(obs, '1236;91347;543;561')
        obs = lineage_str('561', tree, self.proteo['names'])
        exp = 'Gammaproteobacteria;Enterobacterales;Enterobacteriaceae'
        self.assertEqual(obs, exp)
        obs = lineage_str('1224', tree)
        self.assertEqual(obs, '')
        obs = lineage_str('1236', tree)
        self.assertEqual(obs, '')
        obs = lineage_str('0000', tree)
        self.assertEqual(obs, '')

    def test_fill_root(self):
        # already has root
        tree = {'root': 'root',
                'left': 'root',
                'right': 'root'}
        self.assertEqual(fill_root(tree), 'root')

        # multiple crown nodes with parent as self
        tree = {'2': '2',
                '3': '3',
                '4': '2',
                '5': '4',
                '6': '3'}
        self.assertEqual(fill_root(tree), '1')
        self.assertIn('1', tree)
        self.assertEqual(tree['1'], '1')

        # multiple crown nodes with parent as None
        tree = {'1': None,
                '2': '1',
                '3': '1',
                '4': None,
                '5': '4'}
        self.assertEqual(fill_root(tree), '6')
        self.assertIn('6', tree)
        self.assertEqual(tree['1'], '6')
        self.assertEqual(tree['4'], '6')

        # multiple crown nodes with one non-existent parent
        tree = {'Proteobacteria': 'Bacteria',
                'Actinobacteria': 'Bacteria',
                'Firmicutes':     'Bacteria'}
        self.assertEqual(fill_root(tree), 'Bacteria')
        self.assertIn('Bacteria', tree)
        self.assertEqual(tree['Bacteria'], 'Bacteria')

        # multiple crown nodes with multiple non-existent parents
        tree = {'Proteobacteria': 'Bacteria',
                'Actinobacteria': 'Bacteria',
                'Firmicutes':     'Bacteria',
                'Euryarchaeota':  'Archaea',
                'Crenarchaeota':  'Archaea'}
        self.assertEqual(fill_root(tree), '1')
        self.assertIn('1', tree)
        self.assertEqual(tree['Bacteria'], '1')
        self.assertEqual(tree['Archaea'], '1')

    def test_find_rank(self):
        # simple case
        tree = {'1':     '1',
                '2':     '1',
                '1224':  '2',
                '1236':  '1224',
                '28211': '1224',
                '32066': '2'}
        ranks = {'1': 'no rank',
                 '2': 'superkingdom',
                 '1224': 'phylum',
                 '1236': 'class',
                 '28211': 'class',
                 '32066': 'phylum'}
        self.assertEqual(find_rank('1236', 'class', tree, ranks), '1236')
        self.assertEqual(find_rank('1236', 'phylum', tree, ranks), '1224')
        self.assertEqual(find_rank('32066', 'superkingdom', tree, ranks), '2')
        self.assertIsNone(find_rank('12345', 'species', tree, ranks))

        # proteo tree
        args = self.proteo['tree'], self.proteo['ranks']
        # Escherichia coli (species) => Enterobacteriaceae (family)
        self.assertEqual(find_rank('561', 'family', *args), '543')
        # Bordetella (genus) => Betaproteobacteria (class)
        self.assertEqual(find_rank('517', 'class', *args), '28216')
        # Yersinia pestis (species) => self (species)
        self.assertEqual(find_rank('632', 'species', *args), '632')
        # Rickettsiales (order) is already above genus
        self.assertIsNone(find_rank('766', 'genus', *args))
        # no rank is named as "tribe"
        self.assertIsNone(find_rank('317', 'tribe', *args))

    def test_find_lca(self):
        # simple tree
        tsv = ('1	1',
               '2	1',
               '1224	2',
               '1236	1224',
               '28211	1224',
               '32066	2',
               '2157	1')
        tree = read_nodes(tsv)[0]
        self.assertEqual(find_lca(['1236'], tree), '1236')
        self.assertEqual(find_lca(['1236', '32066'], tree), '2')
        self.assertEqual(find_lca(['1236', '2157'], tree), '1')
        self.assertEqual(find_lca(['1'], tree), '1')
        self.assertIsNone(find_lca(['1234'], tree))
        self.assertIsNone(find_lca(['1234', '1'], tree))
        self.assertIsNone(find_lca(['1', '2', '3'], tree))
        self.assertEqual(find_lca(['1', '1'], tree), '1')

        # broken tree, only for unit test coverage
        tree = {'1': '1', '2': '2'}
        self.assertEqual(find_lca(['1', '2'], tree), '1')

        # proteo tree
        # Escherichia, Enterobacter => Enterobacteriaceae
        obs, exp = ['562', '550'], '543'
        self.assertEqual(find_lca(obs, self.proteo['tree']), exp)
        # Helicobacter, Desulfovibrio => Proteobacteria
        obs, exp = ['209', '881'], '1224'
        self.assertEqual(find_lca(obs, self.proteo['tree']), exp)
        # Rickettsia, Sinorhizobium, Rhodobacter => Alphaproteobacteria
        obs, exp = ['782', '382', '1063'], '28211'
        self.assertEqual(find_lca(obs, self.proteo['tree']), exp)


if __name__ == '__main__':
    main()
