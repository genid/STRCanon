#!/usr/bin/env python3
"""Tests for the forensic marker panel and CE allele inference.

Run with:  python -m unittest -v test_markers

The sequences below are real NIST STRSeq records (BioProject PRJNA380127),
each carrying the length-based CE allele that was published with it. They are
the ground truth the bundled panel was fitted against: the anchors and offsets
in str_markers.PANEL_TSV reproduce the published allele for every one of the
1208 STRSeq records they can be placed in, and mis-assign none.
"""

import unittest
from collections import Counter

import strcanon as S
import str_markers as M

# (locus, accession, published CE allele, sequence)
STRSEQ = [
    ('TH01', 'MH085115.2', '6',
     'TGCAGGTCACAGGGAACACAGACTCCATGGTGAATGAATGAATGAATGAATGAATGAGGGAAATAAGGGAG'
     'GAACAGGCCAATGGGAATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCCTATTTGCCCAGGACCC'
     'CCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGCTGCCTCCCCAAGCAGGCAGGCTGGTTGG'
     'GGTGC'),
    # The textbook microvariant: 9 whole AATG repeats plus a 3-base ATG.
    ('TH01', 'MH085124.2', '9.3',
     'TGCAGGTCACAGGGAACACAGACTCCATGGTGAATGAATGAATGAATGAATGAATGATGAATGAATGAATG'
     'AGGGAAATAAGGGAGGAACAGGCCAATGGGAATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCCT'
     'ATTTGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGCTGCCTCCCCAAGC'
     'AGGCAGGCTGGTTGGGGTGC'),
    ('D3S1358', 'OK330002.2', '16',
     'TGCCCACTTCTGCCCAGGGATCTATTTTTCTGTGGTGTGTATTCCCTGTGCCTTTGGGGGCATCTCTTATA'
     'CTCATGAAATCAACAGAGGCTTGCATGTATCTATCTGTCTGTCTCTCTATCTATCTATCTATCTATCTATC'
     'TATCTATCTATCTATCTATCTATGAGACAGGGTCTTGCTCTGTC'),
    ('D13S317', 'MH167239.2', '15',
     'TTCTTTAGTGGGCATCCGTGACTCTCTGGACTCTGACCCATCTAACGCCTATCTGTATTTACAAATACATT'
     'ATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCTATCAATCATCT'
     'ATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATGGCTCAACCCAAGTTGAAGGAGGAGATTT'),
    ('D21S11', 'OK330016.2', '29',
     'AAATATGTGAGTCAATTCCCCAAGTGAATTGCCTTCTATCTATCTATCTATCTATCTATCTGTCTGTCTGT'
     'CTGTCTGTCTGTCTATCTATCTATATCTATCTATCTATCATCTATCTATCCATATCTATCTATCTATCTAT'
     'CTATCTATCTATCTATCTATCGTCTATCTATCCAGTCTATCTACCTCCTATTAGTCT'),
    ('vWA', 'OK330028.2', '17',
     'CATAGGATGGATGGATAGATGGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATAGATGGA'
     'CAGACAGACAGACAGATAGATCAATCCAAGTCACATACTGATTATTCTTATCATCCACTAGGGCTTTCACA'
     'TCTCAGCCAAGTCAACTTGGATCCTCTAGACCTGTTTCTTCTTCTGGAA'),
    ('D22S1045', 'OK623616.1', '18',
     'CGTTGGAATTCCCCAAACTGGCCAGTTCCTCTCCACCCTATAGACCCTGTCCTAGCCTTCTTATAGCTGCT'
     'ATGGGGGCTAGATTTTCCCCGATGATAGTAGTCTCATTATTATTATTATTATTATTATTATTATTATTATT'
     'GTTATTATTACTATTATTGTTATAAAAATATTGCCAAT'),
]

# TH01 MH085118.2, published allele 7, carries a flanking SNP: AACAGAGACT
# where the reference reads AACACAGACT. It is found only because the anchors
# tolerate a mismatch -- with a budget of 0 the locus is missed entirely.
FLANKING_SNP = (
    'MH085118.2', '7',
    'TGCAGGTCACAGGGAACAGAGACTCCATGGTGAATGAATGAATGAATGAATGAATGAATGAGGGAAATAAG'
    'GGAGGAACAGGCCAATGGGAATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCCTATTTGCCCAGG'
    'ACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGCTGCCTCCCCAAGCAGGCAGGCTGG'
    'TTGGGGTGC')

# D21S11 MZ325989.1 is trimmed so tightly that the 5' anchor runs off the start
# of the record. There is no flanking sequence left to anchor to, so there is no
# CE call to be made -- the sequence is still named, just not sized.
NO_FLANK = (
    'MZ325989.1', '25',
    'TGAATTGCCTTCTATCTATCTATCTATCTGTCTGTCTGTCTGTCTGTCTGTCTATCTATCTATATCTATCT'
    'ATCTATCATCTATCTATCCATATCTATCTATCTATCTATCTATCTATCTATCGTCTATCTATCCAGTCTAT'
    'CTACCTCCTATTAGTCT')

# Same locus, same CE allele, different molecule. MH085123.2 opens with a TGAG
# repeat where MH085122.2 has TGAA, so one reads 8 whole AATG units + 3 bases
# and the other 9 + 3. Identical region length, identical CE call, different
# sequence -- exactly the pair CE alone cannot separate.
ISO_ALLELE = [
    ('MH085123.2', '9',
     'TGCAGGTCACAGGGAACACAGACTCCATGGTGAGTGAATGAATGAATGAATGAATGAATGAATGAATGAGG'
     'GAAATAAGGGAGGAACAGGCCAATGGGAATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCCTATT'
     'TGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGCTGCCTCCCCAAGCAGG'
     'CAGGCTGGTTGGGGTGC'),
    ('MH085122.2', '9',
     'TGCAGGTCACAGGGAACACAGACTCCATGGTGAATGAATGAATGAATGAATGAATGAATGAATGAATGAGG'
     'GAAATAAGGGAGGAACAGGCCAATGGGAATCACCCCAGAGCCCAGATACCCTTTGAATTTTGCCCCCTATT'
     'TGCCCAGGACCCCCCACCATGAGCTGCTGCTAGAGCCTGGGAAGGGCCTTGGGGTTGCCTCCCCAAGCAGG'
     'CAGGCTGGTTGGGGTGC'),
]

# A D13S317 duplication allele: the whole repeat block *and its 3' flank* occur
# twice, so the 3' anchor matches inside the array as well as at the end.
DUPLICATION = (
    'MN983127.1', '28.2',
    'TCTGACCCATCTAACGCCTATCTGTATTTACAAATACATTATCTATCTATCTATCTATCTATCTATCTATC'
    'TATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGGCTGCCTATATCTATCTATCTATCT'
    'ATCTATCTATCAATCAATCATCTATCTATCTTTCTGTCTGTCTTTTTGGG')


def nom(seq):
    return S.build_nomenclature(seq, S.call_stretches(seq, 3), False, True, True)[0]


class TestAlleleArithmetic(unittest.TestCase):
    """region_len = offset + period * allele, and back again."""

    def test_whole_and_microvariant_alleles(self):
        th01 = next(m for m in M.DEFAULT_PANEL if m.name == 'TH01')
        # offset 3, period 4: allele 9 spans 3 + 36 = 39 bases, 9.3 spans 42.
        self.assertEqual(str(th01.allele(39)), '9')
        self.assertEqual(str(th01.allele(42)), '9.3')
        self.assertEqual(str(th01.allele(43)), '10')

    def test_allele_and_region_len_are_inverses(self):
        for m in M.DEFAULT_PANEL:
            for repeats in range(0, 40):
                for rem in range(m.period):
                    a = m.allele(m.region_len(repeats, rem))
                    self.assertEqual((a.repeats, a.remainder), (repeats, rem),
                                     f'{m.name} {repeats}.{rem}')

    def test_region_shorter_than_a_zero_allele_is_not_an_allele(self):
        th01 = next(m for m in M.DEFAULT_PANEL if m.name == 'TH01')
        self.assertIsNone(th01.allele(th01.offset - 1))

    def test_allele_formatting(self):
        self.assertEqual(str(M.CEAllele(12, 0)), '12')
        self.assertEqual(str(M.CEAllele(15, 3)), '15.3')


class TestPanel(unittest.TestCase):

    def test_bundled_panel_parses_and_covers_autosomal_x_and_y(self):
        self.assertEqual(len(M.DEFAULT_PANEL), 63)
        names = {m.name for m in M.DEFAULT_PANEL}
        for expected in ('TH01', 'TPOX', 'CSF1PO', 'FGA', 'vWA', 'D21S11',
                         'PentaD', 'PentaE', 'SE33',    # extra autosomal
                         'DXS7132', 'HPRTB',            # X
                         'DYS391', 'DYS438', 'DYS390',  # Y
                         'DYS570', 'DYS626', 'DYS711',  # RM Y-STRs (RMplex)
                         'DYS712', 'DYF403S1b'):
            self.assertIn(expected, names)

        counts = Counter(m.mtype for m in M.DEFAULT_PANEL)
        self.assertEqual(counts, Counter({'AUTOSOMAL': 26, 'X': 5, 'Y': 32}))

    def test_anchors_are_dna_and_periods_sane(self):
        for m in M.DEFAULT_PANEL:
            self.assertTrue(set(m.flank5) <= set('ACGT'), m.name)
            self.assertTrue(set(m.flank3) <= set('ACGT'), m.name)
            self.assertIn(m.period, (3, 4, 5, 6), m.name)
            self.assertGreaterEqual(m.offset, 0, m.name)

    def test_exported_tsv_matches_the_embedded_panel(self):
        with open('str_markers.tsv') as fh:
            exported = M.parse_panel(fh.read())
        self.assertEqual([(m.name, m.flank5, m.flank3, m.period, m.offset)
                          for m in exported],
                         [(m.name, m.flank5, m.flank3, m.period, m.offset)
                          for m in M.DEFAULT_PANEL])

    def test_a_short_row_is_rejected(self):
        with self.assertRaises(ValueError):
            M.parse_panel('X\tAUTOSOMAL\tACGT\tACGT\n')

    def test_comments_and_blank_lines_are_ignored(self):
        panel = M.parse_panel(
            '# a comment\n\nX\tAUTOSOMAL\tAAAA\tCCCC\tAGAT\t4\t0\n')
        self.assertEqual(len(panel), 1)
        self.assertEqual(panel[0].name, 'X')


class TestAnchorMatching(unittest.TestCase):

    def test_exact_match_is_leftmost(self):
        self.assertEqual(M.find_anchor('TTACGTTT', 'ACGT', 0), (2, 0, 1))

    def test_mismatch_budget_is_respected(self):
        self.assertIsNone(M.find_anchor('TTACCTTT', 'ACGT', 0))
        self.assertEqual(M.find_anchor('TTACCTTT', 'ACGT', 1), (2, 1, 1))

    def test_an_exact_match_beats_an_earlier_mismatched_one(self):
        # ACCT (1 mismatch) at 0, ACGT (exact) at 4 -- the exact one wins.
        self.assertEqual(M.find_anchor('ACCTACGT', 'ACGT', 1), (4, 0, 1))

    def test_repeated_anchor_is_reported_as_ambiguous(self):
        self.assertEqual(M.find_anchor('ACGTACGT', 'ACGT', 0), (0, 0, 2))

    def test_iupac_codes_match_their_base_sets(self):
        self.assertEqual(M.find_anchor('TTACGT', 'AYGT', 0), (2, 0, 1))  # Y = C/T
        self.assertIsNone(M.find_anchor('TTAAGT', 'AYGT', 0))


class TestIdentify(unittest.TestCase):

    def test_published_strseq_alleles_are_reproduced(self):
        for locus, acc, allele, seq in STRSEQ:
            with self.subTest(record=acc):
                hit = M.identify(seq)
                self.assertIsNotNone(hit, f'{acc}: no marker recognised')
                self.assertEqual(hit.marker.name, locus, acc)
                self.assertEqual(str(hit.allele), allele, acc)

    def test_the_same_call_on_either_strand(self):
        for locus, acc, allele, seq in STRSEQ:
            with self.subTest(record=acc):
                hit = M.identify(S.reverse_complement(seq))
                self.assertEqual(hit.marker.name, locus, acc)
                self.assertEqual(str(hit.allele), allele, acc)
                self.assertEqual(hit.strand, '-', acc)

    def test_a_flanking_snp_is_absorbed_by_the_mismatch_budget(self):
        acc, allele, seq = FLANKING_SNP
        hit = M.identify(seq, max_mismatches=1)
        self.assertEqual(hit.marker.name, 'TH01', acc)
        self.assertEqual(str(hit.allele), allele, acc)
        self.assertEqual(hit.mismatches, 1, acc)
        # With no budget the SNP hides the anchor and the locus is missed.
        self.assertIsNone(M.identify(seq, max_mismatches=0), acc)

    def test_a_sequence_trimmed_past_its_anchors_gets_no_call(self):
        """No flanking sequence, no CE allele -- the length is not recoverable."""
        acc, _, seq = NO_FLANK
        self.assertIsNone(M.identify(seq), acc)

    def test_a_duplication_allele_is_refused_not_guessed(self):
        """An anchor matching twice leaves no single region -- so, no call.

        Calling this one leftmost would confidently report allele 10 for a
        sequence whose published allele is 28.2.
        """
        acc, allele, seq = DUPLICATION
        hit = M.identify(seq)
        self.assertIsNotNone(hit)
        self.assertEqual(hit.marker.name, 'D13S317')
        self.assertTrue(hit.ambiguous, acc)
        self.assertIsNone(hit.allele, acc)

    def test_sequence_without_anchors_gets_no_marker(self):
        self.assertIsNone(M.identify('AGAT' * 12))

    def test_only_restricts_the_search(self):
        _, _, allele, seq = STRSEQ[0]
        self.assertEqual(str(M.identify(seq, only='TH01').allele), allele)
        self.assertIsNone(M.identify(seq, only='FGA'))
        with self.assertRaises(KeyError):
            M.identify(seq, only='NOSUCHMARKER')


class TestIsoAlleles(unittest.TestCase):
    """One CE allele, two molecules -- the reason the pair is reported."""

    def test_same_ce_allele_but_different_nomenclature(self):
        (acc_a, allele_a, seq_a), (acc_b, allele_b, seq_b) = ISO_ALLELE
        hit_a, hit_b = M.identify(seq_a), M.identify(seq_b)

        # CE cannot separate them: same locus, same allele, same length.
        self.assertEqual(str(hit_a.allele), str(hit_b.allele))
        self.assertEqual(allele_a, allele_b)
        self.assertEqual(hit_a.region_len, hit_b.region_len)

        # The sequence does: different molecules, different brackets.
        self.assertNotEqual(seq_a, seq_b)
        self.assertNotEqual(nom(seq_a), nom(seq_b))


class TestAnnotation(unittest.TestCase):

    def test_annotation_reports_marker_allele_and_lengths(self):
        _, _, _, seq = STRSEQ[1]              # TH01 9.3
        ann = S.ce_annotation(seq, M.identify(seq))
        self.assertEqual(ann, f'{{TH01 CE=9.3 region=42bp len={len(seq)}bp}}')

    def test_unrecognised_sequence_still_reports_its_length(self):
        seq = 'AGAT' * 12
        self.assertEqual(S.ce_annotation(seq, M.identify(seq)),
                         f'{{len={len(seq)}bp}}')

    def test_ambiguous_hit_says_so_instead_of_giving_a_number(self):
        _, _, seq = DUPLICATION
        ann = S.ce_annotation(seq, M.identify(seq))
        self.assertIn('CE=ambiguous', ann)
        self.assertIn('D13S317', ann)

    def test_annotated_line_still_round_trips_through_expand(self):
        """The braces sit outside the repeat grammar, so --expand ignores them."""
        for _, acc, _, seq in STRSEQ:
            with self.subTest(record=acc):
                stretches = S.call_stretches(seq, 3)
                body = S.build_nomenclature(seq, stretches, False, True, False)[0]
                line = f'{S.ce_annotation(seq, M.identify(seq))} {body}'
                self.assertEqual(S.convert_nomenclature_to_sequence(line), seq)


if __name__ == '__main__':
    unittest.main()
