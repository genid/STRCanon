#!/usr/bin/env python3
"""str_markers.py — forensic marker panel and CE allele inference.

A capillary electrophoresis (CE) allele is a *length* call: the fragment is
sized, and its size is read off as a repeat count. Sequence-based nomenclature
carries far more information, but it does not by itself say which CE allele a
sequence would have been called as -- and CE allele numbers are what the rest
of forensic practice (databases, match reports, published frequencies) is
written in. So STRCanon reports both, and the pair makes iso-alleles obvious:
two sequences with the same CE allele but different bracket structures are
different molecules that CE cannot tell apart.

Inferring the CE allele needs a reference panel. Each marker supplies:

    flank5, flank3   anchor sequences bracketing the STR region
    period           repeat unit length (4 = tetranucleotide, 5 = penta, ...)
    offset           number of *non-repeat* bases between the anchors,
                     i.e. the region length of a hypothetical allele 0

The region between the anchors is then

    region_len = offset + period * allele

so a sequence whose inter-anchor region measures `region_len` bases is

    allele, remainder = divmod(region_len - offset, period)

reported as `allele` when the remainder is 0 and `allele.remainder` otherwise
(TH01 9.3 = 9 full AATG repeats + 3 leftover bases).

The anchors are what make this kit-independent. An amplicon's absolute size
depends on where a kit's primers sit, but the offset is calibrated against the
*same* anchor pair used to measure the region, so the primer-dependent constant
cancels: any sequence containing both anchors yields the same CE allele no
matter how much flanking sequence came with it.

Panel format (tab-separated, '#' comments):

    marker  type  flank5  flank3  motif  period  offset

which is the same column layout STRait Razor uses for its locus configs, so a
lab's existing `.config` file can be passed to `--markers` directly.
"""

import os
import sys

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}

# IUPAC ambiguity codes, so a panel may use a degenerate base at a known SNP
# position in an anchor rather than spending its mismatch budget on it.
IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT', 'K': 'GT', 'M': 'AC',
    'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT',
}


def reverse_complement(seq: str) -> str:
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(seq))


# ── The panel ─────────────────────────────────────────────────────────────────
#
# The 20 CODIS core loci. Anchors and offsets are *derived from and validated
# against* NIST STRSeq (BioProject PRJNA380127), the public catalogue of
# sequenced STR alleles: 1208 GenBank records carrying a published, length-based
# CE allele. For every record these anchors can be placed, the arithmetic above
# reproduces that published allele exactly -- see test_markers.py, which asserts
# it, and the README for how the panel was fitted.
#
# The anchors sit as close to the repeat array as they can while staying
# unambiguous. That matters: anchors further out (as used by read-based tools,
# which see whole amplicons) fall outside the ISFG reported range that STRSeq
# and similar catalogues are trimmed to, and simply are not present in the data.
#
# Panels are swappable -- `--markers FILE` takes this same 7-column layout,
# which is also STRait Razor's locus-config format, so a lab's existing kit
# config (ForenSeq, PowerSeq, Y-STR, SE33, ...) can be passed in directly.

PROVENANCE = ('fitted to NIST STRSeq (BioProject PRJNA380127); '
              '26 autosomal, 5 X and 20 Y loci')

PANEL_TSV = """\n#marker	type	flank5	flank3	motif	period	offset
CSF1PO	AUTOSOMAL	GACCCTGTTCTAAGTACTTC	TCTATCTATGAAGGCAGTTA	AGAT	4	15
D10S1248	AUTOSOMAL	TTGAACAAATGAGTGAGT	ATGAAGACAATACAACCAGA	AAGG	4	0
D12ATA63	AUTOSOMAL	GCAATTTAAAAA	CTTGAGACAGGGTCTCGCTC	TAT	3	2
D12S391	AUTOSOMAL	ATCAATGGATGCATAGGTAG	GAGGGGATTTATTAGAGGAA	ACAG	4	0
D13S317	AUTOSOMAL	ATCTGTATTTACAAATACAT	TCTGTCTGTCTTTTTGGG	AGAT	4	21
D16S539	AUTOSOMAL	TACAGACAGACAGACAGGTG	TCATTGAAAGACAAAACAGA	AGAT	4	0
D18S51	AUTOSOMAL	ACAAATTGAGACCTTGTCTC	AAGAGAGAGGAAAGAAAGAG	AAAG	4	1
D19S433	AUTOSOMAL	ATAAAAATCTTCTCTCTT	TTTTCCTTCAACAGAATCTT	AAGG	4	11
D1S1656	AUTOSOMAL	ACAATTAAACACACACACAC	ATCATACAGTTGACCCTTGA	AGAT	4	0
D20S482	AUTOSOMAL	GACACCGAACCAATAA	GATTTATTATAGGAATTGAT	AGAT	4	4
D21S11	AUTOSOMAL	ATTCCCCAAGTGAATTGCCT	GTCTATCTACCTCCTATTAG	ACAG	4	26
D22S1045	AUTOSOMAL	CTTATAGCTGCTATGGGGGC	CTATTATTGTTATAAAAATA	AAT	3	19
D2S1338	AUTOSOMAL	TGCAGGAGGGAAGGAAGGAC	TTCTGTTTCCAAATCCACTG	AAGG	4	0
D2S441	AUTOSOMAL	CCAGGAACTGTGGCTCATCT	ATCTATCTATATCA	AGAT	4	0
D3S1358	AUTOSOMAL	AACAGAGGCTTGCATGTATC	AGACAGGGTCTTGCTC	AGAT	4	0
D5S818	AUTOSOMAL	TCTGTATCCTTATTTATACC	TCAAAATATTACGTAAGGAT	AGAT	4	3
D6S1043	AUTOSOMAL	GCCACTTCCCATAATAAATC	GATCTATCAATCTATTGATC	ATCT	4	2
D7S820	AUTOSOMAL	ATTTAGTGAGATTAAAAAAA	GTTAGTTCGTTCTAAACTAT	AGAT	4	14
D8S1179	AUTOSOMAL	TGTGTACATTCG	TCCCCACAGTGAAAATAATC	AGAT	4	3
FGA	AUTOSOMAL	CAAAAAAGAAAGGAAGAAAG	TAGCTTGTAAATATGC	AAAG	4	0
PentaD	AUTOSOMAL	ATCTCAAGAAAGAAAAAAAA	AAGGGGAAAAAAAGAGAATC	TCTTT	5	8
PentaE	AUTOSOMAL	AGAAAACTCCTTACAA	GAGACTGAGTCTTGCTCAGT	TCTTT	5	3
SE33	AUTOSOMAL	ACTTGCTCTTTCTTTCCTTC	TGACGGAGTTTCACTCTTGT	AAAG	4	112
TH01	AUTOSOMAL	AGGGAACACAGACTCCATGG	GGGAAATAAGGGAGGAACAG	AATG	4	3
TPOX	AUTOSOMAL	GGCACTTAGGGAACCCTCAC	TTTGGGCAAATAAACGCT	AATG	4	2
vWA	AUTOSOMAL	CATAGGATGGATGGATAGAT	AGATCAATCCAAGTCACATA	ACAG	4	0
DXS10074	X	TACACACACAGAGAGAGAGA	AAGAAAGAAAGGAAGAAAAT	AAGA	4	3
DXS7132	X	AACCAATAGGATAGATAGAT	AGATGAGAGGGGATTTATTA	TAGA	4	0
DXS7423	X	CAAATAAATGAATGAGTATG	GGAGGAAATCTGGG	TGGA	4	11
DXS8378	X	AAAAAATAAATAAATAAAAT	TGACCTGCCAGGAGCAGGGG	ATAG	4	0
HPRTB	X	CTTTGTCTCTATCTCTATCT	AAGCAAATTCATGCCCTTCT	TAGA	4	11
DYF387S1	Y	GAAGAAAGAGAAAA	ATAAAAAAAACTGTGGTA	AAGG	4	3
DYS19	Y	GGTTAAGGAGAGTGTCACTA	AAACACTATATATATATAAC	TCTA	4	6
DYS389II	Y	GATAGATTGATAGAGGGAGG	ACAGACAGACACACACATAG	TAGA	4	41
DYS390	Y	GTATACTCAGAAACAAGGAA	TAGATAGAATATATTATGGG	GATA	4	14
DYS391	Y	TATCTGTCTGTCTG	GCCTATCTGCCTGCCTACCT	TCTA	4	3
DYS392	Y	ACCAATCCCATTCCTTAGTA	ATAAATGGTGATACAAGAAA	TAT	3	2
DYS437	Y	ATGCCCATCCGG	ATCATCTGTGAATGACAGGG	TCTA	4	75
DYS438	Y	GTAAACAGTATA	ATTTGAAATGGAGTTTCACT	TTTTC	5	1
DYS439	Y	AAGGTGATAGATATAC	AAGTATAAGTAAAGAGATGA	TATC	4	47
DYS448	Y	AGATAGAGACATGGATAA	GGTAAAGATAGAGATAAA	AGAGAT	6	47
DYS456	Y	TGGGACCTTGTGATAATGTA	TTCCATTAGTTCTGTCCCTC	AGAT	4	1
DYS481	Y	CTAACGCTGTTCAGCATGCT	TTTTGAGTCTTG	CTT	3	1
DYS505	Y	CTTTCTCTGTTCTTTTTCTC	TTTCCCTCCTTCTTTCTCTT	TCCT	4	2
DYS533	Y	CATCTTTCTAGCTAGCTATC	ATCTATCATCTTCTATTGTT	AGAT	4	3
DYS549	Y	TGATAGATGATTAGAAAGAT	AAAAATCTACATAAACAAAA	TATC	4	2
DYS570	Y	TGGCTGTGTCCTCCAAGTTC	TTTTTGTAGATAGG	TTTC	4	2
DYS576	Y	CAAGACCTCATCTCTGAATA	AAGCCAAGACAAATACGCTT	CTTT	4	3
DYS612	Y	TTGCCTCCTCCTCCTCCTCT	TTTTCTTTTGCCTTCCCTCA	TCT	3	0
DYS635	Y	TGAATGGATAAAGAAAATGT	GATTCTATGCAAAGTGAGAA	TAGA	4	2
DYS643	Y	GCCTGGTTAAACTACTGTGC	CTTTTTAAAACTTTTTACTT	CTTTT	5	4
"""


# ── Marker / allele types ─────────────────────────────────────────────────────

class Marker:
    __slots__ = ('name', 'mtype', 'flank5', 'flank3', 'motif', 'period', 'offset')

    def __init__(self, name, mtype, flank5, flank3, motif, period, offset):
        self.name = name
        self.mtype = mtype
        self.flank5 = flank5
        self.flank3 = flank3
        self.motif = motif
        self.period = int(period)
        self.offset = int(offset)
        if self.period < 1:
            raise ValueError(f'{name}: period must be >= 1')
        if self.offset < 0:
            raise ValueError(f'{name}: offset must be >= 0')

    def allele(self, region_len: int):
        """CE allele for an inter-anchor region of `region_len` bases.

        None when the region is shorter than the offset, i.e. shorter than a
        hypothetical allele 0 -- that is not a real allele, it means the anchors
        matched something they should not have.
        """
        span = region_len - self.offset
        if span < 0:
            return None
        repeats, remainder = divmod(span, self.period)
        return CEAllele(repeats, remainder)

    def region_len(self, repeats: int, remainder: int = 0) -> int:
        """Inverse of allele(): the region length a given CE allele implies."""
        return self.offset + self.period * repeats + remainder

    def __repr__(self):
        return f'<Marker {self.name} period={self.period} offset={self.offset}>'


class CEAllele:
    """A CE allele: `repeats` full units plus `remainder` leftover bases."""

    __slots__ = ('repeats', 'remainder')

    def __init__(self, repeats: int, remainder: int):
        self.repeats = repeats
        self.remainder = remainder

    def __str__(self):
        return (f'{self.repeats}.{self.remainder}' if self.remainder
                else str(self.repeats))

    def __eq__(self, other):
        if isinstance(other, CEAllele):
            return (self.repeats, self.remainder) == (other.repeats, other.remainder)
        return str(self) == other

    def __hash__(self):
        return hash((self.repeats, self.remainder))

    def __repr__(self):
        return f'CEAllele({self!s})'


class MarkerHit:
    """A marker recognised in a sequence, and the CE allele it implies.

    `allele` is None when `ambiguous` is set: an anchor matched in more than one
    place, so the STR region has no single length and no call can be made.
    """

    __slots__ = ('marker', 'allele', 'strand', 'region_start', 'region_end',
                 'region', 'mismatches', 'has_motif', 'ambiguous')

    def __init__(self, marker, allele, strand, region_start, region_end, region,
                 mismatches, has_motif, ambiguous=False):
        self.marker = marker
        self.allele = allele
        self.strand = strand              # '+' as given, '-' if found on the RC
        self.region_start = region_start  # inter-anchor region, on `strand`
        self.region_end = region_end
        self.region = region
        self.mismatches = mismatches      # total, across both anchors
        self.has_motif = has_motif        # marker's motif seen in the region
        self.ambiguous = ambiguous        # an anchor matched in several places

    @property
    def region_len(self):
        return len(self.region)

    def __repr__(self):
        return (f'<MarkerHit {self.marker.name} '
                f'{"ambiguous" if self.ambiguous else self.allele} '
                f'{self.region_len}bp {self.strand}>')


# ── Panel parsing ─────────────────────────────────────────────────────────────

def parse_panel(text: str) -> list:
    """Parse a tab-separated panel; ignores blank lines and '#' comments."""
    markers = []
    for lineno, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split('\t')
        if len(fields) < 7:
            raise ValueError(
                f'panel line {lineno}: expected 7 tab-separated fields '
                f'(marker, type, flank5, flank3, motif, period, offset), '
                f'got {len(fields)}')
        try:
            markers.append(Marker(*(f.strip() for f in fields[:7])))
        except ValueError as exc:
            raise ValueError(f'panel line {lineno}: {exc}') from None
    return markers


def load_panel(path=None) -> list:
    """Load the bundled panel, or a user panel / STRait Razor .config."""
    if path is None:
        return parse_panel(PANEL_TSV)
    with open(path) as fh:
        return parse_panel(fh.read())


DEFAULT_PANEL = parse_panel(PANEL_TSV)


# ── Anchor matching ───────────────────────────────────────────────────────────

def _matches_at(seq: str, probe: str, i: int, budget: int):
    """Mismatches of `probe` against seq[i:], or None if over `budget`."""
    mm = 0
    for j, expected in enumerate(probe):
        if seq[i + j] not in IUPAC.get(expected, expected):
            mm += 1
            if mm > budget:
                return None
    return mm


def find_anchor(seq: str, probe: str, max_mismatches: int, start: int = 0):
    """Best match of `probe` in seq[start:]: fewest mismatches, then leftmost.

    Returns (index, mismatches, n_best), where n_best is how many positions tie
    for that best score. n_best > 1 means the anchor lands in more than one
    place and the region it is supposed to delimit is not well defined -- see
    match_marker, which refuses to call rather than picking one.
    """
    if not probe or len(probe) > len(seq) - start:
        return None

    if set(probe) <= {'A', 'C', 'G', 'T'}:
        # Fast path: nothing beats an exact match, so only count those.
        hits = []
        i = seq.find(probe, start)
        while i != -1:
            hits.append(i)
            i = seq.find(probe, i + 1)   # overlapping occurrences count too
        if hits:
            return (hits[0], 0, len(hits))
        if max_mismatches == 0:
            return None

    best_mm = None
    hits = []
    for i in range(start, len(seq) - len(probe) + 1):
        mm = _matches_at(seq, probe, i, max_mismatches)
        if mm is None:
            continue
        if best_mm is None or mm < best_mm:
            best_mm, hits = mm, [i]
        elif mm == best_mm:
            hits.append(i)
    if best_mm is None:
        return None
    return (hits[0], best_mm, len(hits))


def match_marker(seq: str, marker: Marker, max_mismatches: int = 1):
    """Locate `marker` in `seq` (as given, no strand search). None if absent.

    An anchor that matches in several places leaves the STR region undefined,
    so the hit is flagged `ambiguous` and carries no allele. That is not a
    hypothetical: STRSeq MN983127.1 is a D13S317 duplication allele (28.2) in
    which the whole repeat block *and its 3' flank* occur twice, so the 3'
    anchor matches once inside the array. Taking the first match would report a
    confident allele 10 for a sequence that is really a 28.2 -- far worse than
    reporting nothing.
    """
    hit5 = find_anchor(seq, marker.flank5, max_mismatches)
    if hit5 is None:
        return None
    region_start = hit5[0] + len(marker.flank5)

    hit3 = find_anchor(seq, marker.flank3, max_mismatches, start=region_start)
    if hit3 is None:
        return None
    region_end = hit3[0]

    ambiguous = hit5[2] > 1 or hit3[2] > 1
    region = seq[region_start:region_end]
    allele = None if ambiguous else marker.allele(len(region))
    if allele is None and not ambiguous:
        return None

    return MarkerHit(marker, allele, '+', region_start, region_end, region,
                     hit5[1] + hit3[1], marker.motif in region, ambiguous)


def identify(seq: str, panel=None, max_mismatches: int = 1, only: str = None):
    """Best marker hit for `seq`, searching both strands. None if no marker fits.

    Ranked by: a clean call beats an ambiguous one, then fewest anchor
    mismatches, then whether the marker's own motif was seen in the region (a
    hit without it is more likely to be spurious), then name, so the choice is
    deterministic.
    """
    if panel is None:
        panel = DEFAULT_PANEL
    if only is not None:
        wanted = only.casefold()
        panel = [m for m in panel if m.name.casefold() == wanted]
        if not panel:
            raise KeyError(only)

    rc = reverse_complement(seq)
    best = None
    for marker in panel:
        for strand, target in (('+', seq), ('-', rc)):
            hit = match_marker(target, marker, max_mismatches)
            if hit is None:
                continue
            hit.strand = strand
            key = (hit.ambiguous, hit.mismatches, not hit.has_motif, marker.name)
            if best is None or key < best[0]:
                best = (key, hit)
    return best[1] if best else None


# ── TSV export ────────────────────────────────────────────────────────────────

TSV_HEADER = f"""\
# STRCanon forensic marker panel.
#
# region_len = offset + period * allele, measured between flank5 and flank3;
# so  allele, remainder = divmod(region_len - offset, period).
#
# Source: {PROVENANCE}
# Validate against your own kit before casework. Panels are swappable:
# `--markers FILE` reads this same layout. Regenerate with:
#     python str_markers.py --export
#
#marker\ttype\tflank5\tflank3\tmotif\tperiod\toffset
"""


def export_tsv(path='str_markers.tsv') -> int:
    rows = [f'{m.name}\t{m.mtype}\t{m.flank5}\t{m.flank3}\t{m.motif}\t'
            f'{m.period}\t{m.offset}' for m in DEFAULT_PANEL]
    with open(path, 'w', newline='\n') as fh:
        fh.write(TSV_HEADER)
        fh.write('\n'.join(rows) + '\n')
    return len(rows)


def main():
    if '--export' in sys.argv:
        target = 'str_markers.tsv'
        n = export_tsv(target)
        print(f'wrote {n} markers to {target}')
        return

    print(f'{len(DEFAULT_PANEL)} markers ({PROVENANCE})\n')
    print(f"{'marker':<12} {'type':<10} {'period':>6} {'offset':>6}  "
          f"{'allele 10 region':>16}")
    for m in DEFAULT_PANEL:
        print(f'{m.name:<12} {m.mtype:<10} {m.period:>6} {m.offset:>6}  '
              f'{m.region_len(10):>13} bp')


if __name__ == '__main__':
    main()
