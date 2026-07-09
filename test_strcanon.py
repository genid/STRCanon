#!/usr/bin/env python3
"""Tests for strcanon. Run with:  python -m unittest -v test_strcanon"""

import random
import unittest

import strcanon as S


def nom(seq, min_repeats=3, hide_n=False, full_seq_gaps=True, strip_ends_n=False):
    """Nomenclature under lossless settings (flanks kept, gaps written out)."""
    stretches = S.call_stretches(seq, min_repeats)
    return S.build_nomenclature(seq, stretches, hide_n, full_seq_gaps, strip_ends_n)[0]


def rendered_blocks(seq, min_repeats=3):
    """(canonical, repeats, partial) for the blocks build_nomenclature emits.

    find_stretches also returns nested candidates that the cursor in
    build_nomenclature skips; those never reach the user.
    """
    emitted, cursor = [], 0
    for s in S.call_stretches(seq, min_repeats):
        if s['end'] <= cursor:
            continue
        emitted.append((s['canonical'], s['repeats'], s['partial_len']))
        cursor = s['end'] + s['partial_len']
    return emitted


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
        for strand in (seq, S.reverse_complement(seq)):
            rendered = nom(strand)
            self.assertEqual(S.convert_nomenclature_to_sequence(rendered), strand,
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
    """Stretches are called on one strand of the duplex, so the two strands of a
    locus get the same blocks in mirror-image order -- never different calls."""

    def assert_mirrored(self, seq, label=''):
        fwd = rendered_blocks(seq)
        rev = rendered_blocks(S.reverse_complement(seq))
        self.assertEqual(fwd, list(reversed(rev)),
                         f'{label or seq}\n  + {nom(seq)}\n  - {nom(S.reverse_complement(seq))}')

    def test_th01_9_3_names_the_same_blocks_on_either_strand(self):
        # The 3-base ATG interruption used to be absorbed as 6.1/3.2 on the plus
        # strand but 3.3/6 on the minus strand.
        seq = LOCI['TH01 a9.3']
        self.assertEqual(rendered_blocks(seq), [('AATG', 6, 1), ('AATG', 3, 2)])
        self.assertEqual(rendered_blocks(S.reverse_complement(seq)),
                         [('AATG', 3, 2), ('AATG', 6, 1)])

    def test_forensic_loci(self):
        for label, seq in LOCI.items():
            with self.subTest(locus=label):
                self.assert_mirrored(seq, label)

    def test_abutting_runs_of_the_same_motif_class(self):
        """The case that breaks a left-to-right scan: two phases of one motif."""
        rng = random.Random(20240105)
        for _ in range(150):
            k = rng.randint(2, 4)
            m = ''.join(rng.choice('ACGT') for _ in range(k))
            seq = m * rng.randint(3, 6) + (m[1:] + m[0]) * rng.randint(3, 6)
            self.assert_mirrored(seq)

    def test_interrupted_runs(self):
        rng = random.Random(20240106)
        for _ in range(150):
            m = ''.join(rng.choice('ACGT') for _ in range(rng.randint(2, 6)))
            gap = ''.join(rng.choice('ACGT') for _ in range(rng.randint(1, 3)))
            self.assert_mirrored(m * rng.randint(3, 6) + gap + m * rng.randint(3, 6))

    def test_random_sequences(self):
        rng = random.Random(20240103)
        for _ in range(200):
            self.assert_mirrored(''.join(rng.choice('ACGT')
                                         for _ in range(rng.randint(10, 60))))

    def test_palindromic_sequence_is_its_own_reverse_complement(self):
        for core in ('AATT', 'ACGT', 'AGCT'):
            seq = core * 6
            self.assertEqual(seq, S.reverse_complement(seq))
            self.assertEqual(nom(seq), f'[{core}]6')

    def test_analysis_strand_is_the_same_for_both_strands(self):
        rng = random.Random(20240107)
        for _ in range(100):
            seq = ''.join(rng.choice('ACGT') for _ in range(rng.randint(6, 40)))
            self.assertEqual(S.canonical_strand(seq)[0],
                             S.canonical_strand(S.reverse_complement(seq))[0])

    def test_mirror_stretch_preserves_repeats_and_partial(self):
        """rc(w*r + w[:p]) == u*r + u[:p] with u = rot(rc(w), k-p)."""
        rng = random.Random(20240108)
        for seq in list(LOCI.values()) + [
                ''.join(rng.choice('ACGT') for _ in range(rng.randint(10, 50)))
                for _ in range(100)]:
            n = len(seq)
            for s in S.call_stretches(seq, 3):
                m = S.mirror_stretch(s, n)
                self.assertEqual((m['canonical'], m['repeats'], m['partial_len']),
                                 (s['canonical'], s['repeats'], s['partial_len']))
                region = S.reverse_complement(seq)[m['start']:
                                                   m['start'] + m['repeats'] * m['k']
                                                   + m['partial_len']]
                self.assertEqual(region, m['variant'] * m['repeats']
                                 + m['variant'][:m['partial_len']])

    def test_canonical_motifs_are_strand_independent(self):
        """The raw finder already agrees on *which* motifs occur, on either strand."""
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

    def test_reverse_complement_handles_n(self):
        self.assertEqual(S.reverse_complement('AAGNNAAG'), 'CTTNNCTT')

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
