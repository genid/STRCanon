#!/usr/bin/env python3
"""Tests for strcanon. Run with:  python -m unittest -v test_strcanon"""

import random
import unittest

import strcanon as S


def nom(seq, min_repeats=3, hide_n=False, full_seq_gaps=True, strip_ends_n=False):
    """Nomenclature under lossless settings (flanks kept, gaps written out)."""
    stretches = S.find_stretches(seq, min_repeats)
    return S.build_nomenclature(seq, stretches, hide_n, full_seq_gaps, strip_ends_n)[0]


# Sequence-based renderings of common forensic loci, built from their published
# repeat structure. These exercise compound motifs, interruptions and overlaps.
LOCI = {
    'D3S1358 a17': 'TCTA' + 'TCTG' * 2 + 'TCTA' * 14,
    'vWA a17': 'TCTA' + 'TCTG' * 4 + 'TCTA' * 11 + 'TCCA' + 'TCTA',
    'D12S391 a18': 'AGAT' * 11 + 'AGAC' * 6 + 'AGAT',
    'FGA a22': 'TTTC' * 3 + 'TTTT' + 'TTCT' + 'CTTT' * 14 + 'CTCC' + 'TTCC' * 2,
    'D8S1179 a13': 'TCTA' * 13,
    'TH01 a9.3': 'AATG' * 6 + 'ATG' + 'AATG' * 3,
    'D21S11 core': ('TCTA' * 4 + 'TCTG' * 6 + 'TCTA' * 3 + 'TA' + 'TCTA' * 3
                    + 'TCA' + 'TCTA' * 2 + 'TCCA' + 'TA' + 'TCTA' * 11),
}


class TestCanonicalisation(unittest.TestCase):

    def test_matches_the_generated_lookup_table(self):
        """The in-memory lookup must agree with str_canonical_motifs.tsv."""
        with open('str_canonical_motifs.tsv') as fh:
            next(fh)
            rows = [line.rstrip('\n').split('\t') for line in fh if line.strip()]
        self.assertTrue(rows, 'lookup table is empty')

        for motif, canonical, is_rc, k in rows:
            self.assertIn(motif, S.MOTIF_TO_CANONICAL, motif)
            self.assertEqual(S.MOTIF_TO_CANONICAL[motif],
                             (canonical, is_rc == '1'), motif)
            self.assertEqual(int(k), len(motif))
        self.assertEqual(len(rows), len(S.MOTIF_TO_CANONICAL))

    def test_canonical_is_lexicographic_minimum_of_its_class(self):
        for canonical, variants in S.CANONICAL_TO_VARIANTS.items():
            self.assertEqual(canonical, min(variants))

    def test_every_variant_maps_back_to_its_canonical(self):
        for canonical, variants in S.CANONICAL_TO_VARIANTS.items():
            for v in variants:
                self.assertEqual(S.MOTIF_TO_CANONICAL[v][0], canonical)

    def test_sub_period_motifs_are_excluded(self):
        for degenerate in ('AA', 'ATAT', 'AAA', 'ACGACG', 'AGAG'):
            self.assertNotIn(degenerate, S.CANONICAL_TO_VARIANTS)

    def test_reverse_complement_motifs_share_a_canonical(self):
        # The four rotations of AGAT and of its RC ATCT are one class.
        for motif in ('AGAT', 'GATA', 'ATAG', 'TAGA',
                      'ATCT', 'TCTA', 'CTAT', 'TATC'):
            self.assertEqual(S.MOTIF_TO_CANONICAL[motif][0], 'AGAT')


class TestExpand(unittest.TestCase):
    """`[N]-n` is an overlap marker: it rewinds the cursor n bases."""

    def test_overlap_rewinds_rather_than_duplicating(self):
        # [AAG]3.1 writes AAGAAGAAGA, then [N]-3 backs up onto the trailing AGA.
        self.assertEqual(
            S.convert_nomenclature_to_sequence('[AAG]3.1[N]-3[AG]4'),
            'AAGAAGAAGAGAGAG')

    def test_gap_marker_expands_to_ns_preserving_length(self):
        self.assertEqual(S.convert_nomenclature_to_sequence('[N]3[AAG]3'),
                         'NNNAAGAAGAAG')

    def test_literal_gap_is_emitted_verbatim(self):
        self.assertEqual(S.convert_nomenclature_to_sequence('(TGC)[AAG]4[AG]3.1'),
                         'TGC' + 'AAG' * 4 + 'AG' * 3 + 'A')

    def test_rotation_and_reverse_complement(self):
        self.assertEqual(S.convert_nomenclature_to_sequence('[AAG]3'), 'AAGAAGAAG')
        self.assertEqual(S.convert_nomenclature_to_sequence('[AAG>1]3'), 'AGAAGAAGA')
        self.assertEqual(S.convert_nomenclature_to_sequence('[~AAG]3'), 'CTTCTTCTT')
        self.assertEqual(S.convert_nomenclature_to_sequence('[~AAG>2]4.2'),
                         'TTCTTCTTCTTCTT')

    def test_overlap_past_the_start_is_an_error(self):
        with self.assertRaises(ValueError):
            S.convert_nomenclature_to_sequence('[AG]3[N]-9[AG]3')


class TestRoundTrip(unittest.TestCase):
    """expand(nomenclature(seq)) == seq, under lossless settings."""

    def assert_round_trips(self, seq, label=''):
        rendered = nom(seq)
        self.assertEqual(S.convert_nomenclature_to_sequence(rendered), seq,
                         f'{label or seq}: {rendered}')

    def test_readme_example(self):
        self.assert_round_trips('AAGAAGAAGAGAGAG')

    def test_forensic_loci(self):
        for label, seq in LOCI.items():
            with self.subTest(locus=label):
                self.assert_round_trips(seq, label)

    def test_no_repeats_at_all(self):
        self.assert_round_trips('ACGTACGGTCAGGCTA')

    def test_random_sequences(self):
        rng = random.Random(20240101)
        for _ in range(300):
            seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(6, 80)))
            self.assert_round_trips(seq)

    def test_repeat_dense_sequences(self):
        """Concatenated tandem arrays -- where overlapping stretches arise."""
        rng = random.Random(20240102)
        for _ in range(300):
            seq = ''
            for _ in range(rng.randint(1, 4)):
                motif = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 6)))
                seq += motif * rng.randint(3, 8)
            self.assert_round_trips(seq)


class TestFindStretches(unittest.TestCase):

    def test_mononucleotide_needs_five_repeats(self):
        self.assertEqual(nom('CGCGAAAACGCG'), '(CGCGAAAACGCG)')
        self.assertIn('[A]5', nom('CGCGAAAAACGCG'))

    def test_partial_repeat_is_reported(self):
        self.assertEqual(nom('AGAGAGA'), '[AG]3.1')

    def test_reverse_complement_stretch(self):
        self.assertEqual(nom('TTCTTCTTCTTCTT'), '[~AAG>2]4.2')

    def test_min_repeats_is_respected(self):
        self.assertEqual(nom('AGAGAG', min_repeats=3), '[AG]3')
        self.assertEqual(nom('AGAGAG', min_repeats=4), '(AGAGAG)')

    def test_stretches_are_ordered_and_within_bounds(self):
        for seq in LOCI.values():
            stretches = S.find_stretches(seq, 3)
            for s in stretches:
                self.assertGreaterEqual(s['start'], 0)
                self.assertLessEqual(s['end'] + s['partial_len'], len(seq))
                self.assertEqual((s['end'] - s['start']) % s['k'], 0)
            starts = [s['start'] for s in stretches]
            self.assertEqual(starts, sorted(starts))

    def test_stretch_bases_match_the_underlying_sequence(self):
        for seq in LOCI.values():
            for s in S.find_stretches(seq, 3):
                expected = s['variant'] * s['repeats'] + s['variant'][:s['partial_len']]
                actual = seq[s['start']:s['end'] + s['partial_len']]
                self.assertEqual(actual, expected)


class TestStrandSymmetry(unittest.TestCase):

    def test_canonical_motifs_are_strand_independent(self):
        """Reverse-complementing a sequence must not change which motifs are called."""
        rng = random.Random(20240103)
        for _ in range(200):
            seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(10, 60)))
            fwd = sorted(s['canonical'] for s in S.find_stretches(seq, 3))
            rev = sorted(s['canonical']
                         for s in S.find_stretches(S.reverse_complement(seq), 3))
            self.assertEqual(fwd, rev, seq)


class TestSequenceHelpers(unittest.TestCase):

    def test_prepare_sequence_uppercases_and_drops_invalid_chars(self):
        self.assertEqual(S.prepare_sequence('aagAAG'), ('AAGAAG', False))
        self.assertEqual(S.prepare_sequence('AAG-AAG'), ('AAGAAG', True))

    def test_prepare_sequence_expands_nomenclature(self):
        self.assertEqual(S.prepare_sequence('[AAG]3')[0], 'AAGAAGAAG')

    def test_reverse_complement_is_an_involution(self):
        rng = random.Random(20240104)
        for _ in range(100):
            seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 40)))
            self.assertEqual(S.reverse_complement(S.reverse_complement(seq)), seq)

    def test_has_sub_period(self):
        self.assertTrue(S.has_sub_period('ATAT'))
        self.assertTrue(S.has_sub_period('AAA'))
        self.assertFalse(S.has_sub_period('AGAT'))
        self.assertFalse(S.has_sub_period('A'))


if __name__ == '__main__':
    unittest.main()
