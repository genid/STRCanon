#!/usr/bin/env python3
"""strnom.py — canonical STR nomenclature from the command line.

Detects short tandem repeat (STR) stretches in DNA sequences and renders them
in a canonical bracket nomenclature, e.g.  (TGC)[AAG]4[AG]3.1

The canonical motif for every k-mer (k = 1..6) is derived on the fly using the
lexicographic-minimum rule over each motif's equivalence class (all cyclic
rotations of the motif and of its reverse complement). This is the same
canonicalisation used by the companion generate_canonicals.py generator, so no
external lookup table is needed.

Notation
    [AAG]        canonical motif, no rotation
    [AAG>1]      canonical rotated 1 position left (AAG -> AGA)
    [~AAG]       reverse complement of canonical (AAG -> CTT)
    [~AAG>2]     RC of canonical rotated 2 left  (AAG -> GAA -> TTC)
    ]3           3 complete repeats
    ]3.1         3 complete repeats + 1 partial base
    [N]5         5-base gap between stretches
    [N]-2        2-base overlap between stretches
    (ACGT)       literal gap sequence (with --full-seq-gaps)

A stretch may claim a partial repeat whose bases the next stretch also covers;
[N]-2 records that the following stretch starts 2 bases before the previous one
ended. Expansion therefore writes terms through a cursor: [N]-2 rewinds it two
bases and the next stretch overwrites them.

Counts use plain digits, so any output line can be fed straight back in with
--expand. Reconstruction is exact when nothing was suppressed on the way out,
i.e. with --no-strip-ends-n and the default --full-seq-gaps.

Stretches are always called on one strand of the duplex (see canonical_strand),
so a locus and its reverse complement yield the same motifs, repeat counts and
partial lengths. --orientation only chooses which strand they are displayed on.

Examples
    strnom.py AAGAAGAAGAGAGAG
    strnom.py TGCAAGAAGAAG --no-strip-ends-n          # keep the leading (TGC)
    strnom.py --input seqs.txt --visual --matrix
    echo TTCTTCTTCTTCTT | strnom.py                   # -> [~AAG>2]4.2
    strnom.py --expand '(TGC)[AAG]4[AG]3.1'           # nomenclature -> sequence
"""

import sys
import re
import argparse
from itertools import product

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
ACCEPTED = set('ACGTN')
# Minimum number of repeats to call a stretch, per motif length k.
# Mono-nucleotide runs are only interesting when long, hence 5.
MIN_LEN = {1: 5, 2: 2, 3: 2, 4: 2, 5: 2, 6: 2}


# ── Sequence helpers ──────────────────────────────────────────────────────────

def reverse_complement(seq: str) -> str:
    return ''.join(COMPLEMENT[b] for b in reversed(seq))


def get_rotations(seq: str) -> list:
    return [seq[i:] + seq[:i] for i in range(len(seq))]


def has_sub_period(seq: str) -> bool:
    """True if seq is a tandem repeat of a shorter unit (e.g. ATAT = AT*2)."""
    n = len(seq)
    for period in range(1, n):
        if n % period == 0 and seq == seq[:period] * (n // period):
            return True
    return False


def rotation_offset(variant: str, ref: str) -> int:
    """Return r such that variant == ref[r:] + ref[:r] (0 when equal)."""
    if variant == ref:
        return 0
    return (ref + ref).find(variant)


# ── Canonical lookup (generated in memory) ────────────────────────────────────

def build_lookup(k_max: int = 6):
    """Return (motif_to_canonical, canonical_to_variants).

    motif_to_canonical : motif      -> (canonical, is_rc)
    canonical_to_variants : canonical -> [all representations]
    Sub-period motifs (ATAT, etc.) are excluded, as in the generator.
    """
    motif_to_canonical = {}
    canonical_to_variants = {}

    for k in range(1, k_max + 1):
        processed = set()
        for kmer in sorted(''.join(p) for p in product('ACGT', repeat=k)):
            if kmer in processed or has_sub_period(kmer):
                if has_sub_period(kmer):
                    processed.update(get_rotations(kmer))
                    processed.update(get_rotations(reverse_complement(kmer)))
                continue

            fwd = set(get_rotations(kmer))
            group = fwd | set(get_rotations(reverse_complement(kmer)))
            canonical = min(group)

            for motif in group:
                motif_to_canonical[motif] = (canonical, motif not in fwd)
            canonical_to_variants[canonical] = sorted(group)
            processed.update(group)

    return motif_to_canonical, canonical_to_variants


MOTIF_TO_CANONICAL, CANONICAL_TO_VARIANTS = build_lookup()


# ── Core algorithm: find repeat stretches ─────────────────────────────────────

def find_stretches(seq: str, min_repeats: int) -> list:
    """Find STR stretches, preferring the rotation with the most full repeats.

    Suffix partials are recorded in `partial_len`. Where two stretches of the
    *same* canonical overlap, only the first (best) is kept.
    """
    seen = {}  # (raw_start, raw_end) -> best entry

    for canonical, variants in CANONICAL_TO_VARIANTS.items():
        k = len(canonical)
        min_count = max(MIN_LEN[k], min_repeats)

        for variant in variants:
            if variant * min_count not in seq:
                continue
            is_rc = MOTIF_TO_CANONICAL[variant][1]

            for m in re.finditer(f'(?:{variant}){{{min_count},}}', seq):
                raw_start, raw_end = m.start(), m.end()
                repeats = (raw_end - raw_start) // k

                # Longest leading fragment of `variant` after the last full repeat.
                partial_len = 0
                for plen in range(k - 1, 0, -1):
                    if (raw_end + plen <= len(seq)
                            and seq[raw_end:raw_end + plen] == variant[:plen]):
                        partial_len = plen
                        break

                fwd = reverse_complement(variant) if is_rc else variant
                r = rotation_offset(fwd, canonical)

                key = (raw_start, raw_end)
                coverage = repeats * k + partial_len
                old = seen.get(key)
                if (old is None or repeats > old['repeats']
                        or (repeats == old['repeats']
                            and coverage > old['repeats'] * old['k'] + old['partial_len'])):
                    seen[key] = {
                        'start': raw_start, 'end': raw_end, 'variant': variant,
                        'canonical': canonical, 'is_rc': is_rc, 'r': r, 'k': k,
                        'repeats': repeats, 'partial_len': partial_len,
                    }

    all_stretches = sorted(
        seen.values(),
        key=lambda s: (s['start'], -s['repeats'],
                       -(s['repeats'] * s['k'] + s['partial_len'])),
    )

    filtered = []
    for s in all_stretches:
        if any(f['canonical'] == s['canonical']
               and f['start'] < s['end'] and s['start'] < f['end']
               for f in filtered):
            continue
        filtered.append(s)
    return filtered


# ── Analysis strand ───────────────────────────────────────────────────────────
#
# find_stretches() reads 5'->3': the regex takes the leftmost match and a
# partial repeat is the leftover on a run's 3' end. Reverse-complementing a
# sequence moves that leftover to the 5' end, so calling stretches directly on
# whichever strand happened to be sequenced can place partials differently.
# TH01 allele 9.3 is the textbook case: its 3-base ATG interruption is absorbed
# as [AATG]6.1...3.2 on one strand and [~AATG]3.3...6 on the other.
#
# So stretches are always called on one strand of the duplex, chosen from the
# molecule itself rather than from the input: the lexicographically smaller of
# the sequence and its reverse complement. Both strands of a locus canonicalise
# to the same string, hence to the same blocks. The results are then mirrored
# onto whichever strand is being displayed, which changes only the order of the
# blocks and how each is written (~ and >r) -- never the canonical motif, the
# repeat count or the partial length.

def canonical_strand(seq: str):
    """Return (analysis_seq, flipped) for the strand stretches are called on."""
    rev = reverse_complement(seq)
    return (seq, False) if seq <= rev else (rev, True)


def mirror_stretch(s: dict, n: int) -> dict:
    """Map a stretch onto the opposite strand of a length-n sequence.

    A run reads `w * r + w[:p]`. Its reverse complement is `u * r + u[:p]`
    where u = rc(w) rotated left by k - p, so the repeat count and the partial
    length survive the flip; only the motif's rotation changes.
    """
    k, r, p = s['k'], s['repeats'], s['partial_len']
    rc_w = reverse_complement(s['variant'])
    u = rc_w[k - p:] + rc_w[:k - p]

    start = n - (s['start'] + r * k + p)
    canonical, is_rc = MOTIF_TO_CANONICAL[u]
    fwd = reverse_complement(u) if is_rc else u
    return {'start': start, 'end': start + r * k, 'variant': u,
            'canonical': canonical, 'is_rc': is_rc,
            'r': rotation_offset(fwd, canonical), 'k': k,
            'repeats': r, 'partial_len': p}


def mirror_stretches(stretches: list, n: int) -> list:
    return sorted((mirror_stretch(s, n) for s in stretches),
                  key=lambda s: (s['start'], -s['repeats'],
                                 -(s['repeats'] * s['k'] + s['partial_len'])))


def call_stretches(seq: str, min_repeats: int) -> list:
    """Find stretches on the canonical strand, reported in seq's orientation."""
    analysis_seq, flipped = canonical_strand(seq)
    stretches = find_stretches(analysis_seq, min_repeats)
    return mirror_stretches(stretches, len(seq)) if flipped else stretches


# ── Nomenclature rendering ────────────────────────────────────────────────────

def _stretch_label(s: dict) -> str:
    base = ('~' + s['canonical']) if s['is_rc'] else s['canonical']
    rotation = f">{s['r']}" if s['r'] > 0 else ''
    count = f"{s['repeats']}.{s['partial_len']}" if s['partial_len'] > 0 else str(s['repeats'])
    return f"[{base}{rotation}]{count}"


def build_nomenclature(seq, stretches, hide_n, full_seq_gaps, strip_ends_n):
    """Return (nomenclature_string, segments).

    segments: list of dicts {start, seq, label, canonical} for the visual.
    """
    cursor = 0
    parts = []
    segments = []
    is_first = True

    for s in stretches:
        if s['end'] <= cursor:
            continue

        gap = s['start'] - cursor
        if gap > 0:
            if not hide_n and not (strip_ends_n and is_first):
                gap_seq = seq[cursor:s['start']]
                parts.append(f"({gap_seq})" if full_seq_gaps else f"[N]{gap}")
        elif gap < 0:
            parts.append(f"[N]{gap}")  # overlap: gap is negative, e.g. [N]-2

        is_first = False
        vis_seq = s['variant'] * s['repeats'] + s['variant'][:s['partial_len']]
        parts.append(_stretch_label(s))

        base = ('~' + s['canonical']) if s['is_rc'] else s['canonical']
        rotation = f">{s['r']}" if s['r'] > 0 else ''
        count = f"{s['repeats']}.{s['partial_len']}" if s['partial_len'] > 0 else str(s['repeats'])
        segments.append({'start': s['start'], 'seq': vis_seq,
                         'label': f"[{base}{rotation}]x{count}",
                         'canonical': s['canonical']})
        cursor = s['end'] + s['partial_len']

    if cursor < len(seq) and not hide_n and not strip_ends_n:
        tail = seq[cursor:]
        parts.append(f"({tail})" if full_seq_gaps else f"[N]{len(tail)}")

    return ''.join(parts), segments


TOKEN = re.compile(
    r'\(([ACGTN]+)\)|\[(~?)([ACGTN]+)(?:>(\d+))?\](-?\d+)(?:\.(\d+))?'
)


def convert_nomenclature_to_sequence(nomenclature: str) -> str:
    """Expand a bracket nomenclature back to a plain DNA sequence.

    Terms are written through a cursor rather than concatenated, because
    `[N]-n` is an *overlap*: it rewinds the cursor n bases so the following
    stretch re-states bases the previous one already covered. Concatenating
    would emit those bases twice.

    `[N]n` (a gap with `--no-full-seq-gaps`) expands to n literal Ns: the
    bases are unrecoverable, but the length is preserved.
    """
    out = []   # bases written so far, addressed by `cursor`
    cursor = 0

    for lit, tilde, canonical, rot, count, partial in TOKEN.findall(nomenclature):
        if lit:
            piece = lit
        elif canonical == 'N':
            n = int(count)
            if n < 0:                      # overlap: rewind, then overwrite
                cursor += n
                if cursor < 0:
                    raise ValueError(
                        f'overlap [N]{n} rewinds past the start of the sequence')
                continue
            piece = 'N' * n                # gap of unknown bases
        else:
            r = int(rot) if rot else 0
            p = int(partial) if partial else 0
            rotated = canonical[r:] + canonical[:r]
            variant = reverse_complement(rotated) if tilde == '~' else rotated
            piece = variant * int(count) + variant[:p]

        for base in piece:
            if cursor < len(out):
                out[cursor] = base
            else:
                out.append(base)
            cursor += 1

    return ''.join(out)


# ── ANSI colour ───────────────────────────────────────────────────────────────

class Palette:
    CODES = [203, 214, 220, 78, 45, 33, 141, 205, 118, 208, 39, 170, 76, 51, 213, 226]

    def __init__(self, enabled: bool):
        self.enabled = enabled

    def paint(self, text: str, idx: int, underline: bool = False) -> str:
        if not self.enabled:
            return text
        code = self.CODES[idx % len(self.CODES)]
        u = '4;' if underline else ''
        return f"\x1b[{u}38;5;{code}m{text}\x1b[0m"


# ── Visual alignment block ────────────────────────────────────────────────────

def build_visual(seq, segments, color_index, palette):
    """Monospace alignment: full sequence, then one indented line per stretch."""
    if not segments:
        return seq

    col = max(max(s['start'] + len(s['seq']) for s in segments), len(seq)) + 2

    # Per-base colour index for the top line (None if uncovered, -1 if overlap).
    cover = [None] * len(seq)
    for seg in segments:
        idx = color_index[seg['canonical']]
        for i in range(seg['start'], min(seg['start'] + len(seg['seq']), len(seq))):
            cover[i] = idx if cover[i] is None else -1

    # Top line: runs of identical colour.
    top = []
    i = 0
    while i < len(seq):
        c = cover[i]
        j = i + 1
        while j < len(seq) and cover[j] == c:
            j += 1
        chunk = seq[i:j]
        if c is None:
            top.append(chunk)
        elif c == -1:
            top.append(palette.paint(chunk, 0, underline=True) if palette.enabled else chunk)
        else:
            top.append(palette.paint(chunk, c))
        i = j
    lines = [''.join(top)]

    # One line per stretch, aligned under its position.
    for seg in segments:
        idx = color_index[seg['canonical']]
        indent = ' ' * seg['start']
        pad = ' ' * (col - seg['start'] - len(seg['seq']))
        lines.append(indent + palette.paint(seg['seq'], idx) + pad
                     + palette.paint(seg['label'], idx))
    return '\n'.join(lines)


# ── Pairwise Hamming distance matrix ──────────────────────────────────────────

def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def min_hamming_between(c1: str, c2: str):
    v1 = CANONICAL_TO_VARIANTS.get(c1, [c1])
    v2 = CANONICAL_TO_VARIANTS.get(c2, [c2])
    best = min((hamming(a, b) for a in v1 for b in v2 if len(a) == len(b)),
               default=None)
    return best


def build_matrices(canonicals, palette) -> str:
    """Text tables of minimum Hamming distance between equivalence classes."""
    canonicals = sorted(canonicals)
    if len(canonicals) < 2:
        return ''

    by_k = {}
    for c in canonicals:
        by_k.setdefault(len(c), []).append(c)

    out = ['Pairwise motif distance matrix',
           '(minimum Hamming distance between equivalence classes; grouped by k)']

    for k in sorted(by_k):
        cs = sorted(by_k[k])
        if len(cs) < 2:
            out.append(f"\nk = {k}: only one motif ({cs[0]}) — no comparison.")
            continue

        w = max(k, 2)
        out.append(f"\nk = {k}")
        header = ' ' * (w + 2) + '  '.join(c.rjust(w) for c in cs)
        out.append(header)
        for row in cs:
            cells = []
            for col in cs:
                if row == col:
                    cells.append('-'.rjust(w))
                else:
                    d = min_hamming_between(row, col)
                    txt = ('?' if d is None else str(d)).rjust(w)
                    # Highlight one-mutation neighbours (forensically notable).
                    if palette.enabled and d == 1:
                        txt = palette.paint(txt, 0)
                    cells.append(txt)
            out.append(row.rjust(w) + '  ' + '  '.join(cells))
    return '\n'.join(out)


# ── Input handling ────────────────────────────────────────────────────────────

def prepare_sequence(line: str):
    """Uppercase, expand nomenclature if present, strip invalid chars.

    Returns (sequence, had_invalid_chars).
    """
    s = line.strip().upper()
    if '[' in s or ']' in s or '(' in s:
        s = convert_nomenclature_to_sequence(s)
    cleaned = ''.join(ch for ch in s if ch in ACCEPTED)
    return cleaned, len(cleaned) != len(s)


def gather_inputs(args) -> list:
    lines = list(args.sequences)
    if args.input:
        with open(args.input) as fh:
            lines.extend(fh.read().splitlines())
    if not lines and not sys.stdin.isatty():
        lines.extend(sys.stdin.read().splitlines())
    return [ln for ln in lines if ln.strip()]


# ── main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    p = argparse.ArgumentParser(
        description='Render DNA sequences in canonical STR bracket nomenclature.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split('Examples')[-1] if 'Examples' in __doc__ else None,
    )
    p.add_argument('sequences', nargs='*',
                   help='DNA sequence(s) or nomenclature string(s) to analyse')
    p.add_argument('-i', '--input', metavar='FILE',
                   help='read one sequence/nomenclature per line from FILE')
    p.add_argument('--min-repeats', type=int, default=3,
                   help='minimum repeats to call a stretch (mono-nucleotides '
                        'always require 5; default: 3)')
    p.add_argument('--hide-n', action=argparse.BooleanOptionalAction, default=False,
                   help='suppress all [N] gap/overlap markers')
    p.add_argument('--full-seq-gaps', action=argparse.BooleanOptionalAction, default=True,
                   help='show literal bases in (parentheses) instead of [N]n '
                        '(default: on)')
    p.add_argument('--strip-ends-n', action=argparse.BooleanOptionalAction, default=True,
                   help='suppress only the leading/trailing gap markers '
                        '(default: on)')
    p.add_argument('--orientation', choices=('as-given', 'reverse-complement'),
                   default='as-given',
                   help='strand to report on. Stretches are always called on the '
                        'same strand of the duplex, so this changes only the order '
                        'the motifs are shown in and how each is written, never how '
                        'they are called (default: as-given)')

    adv = p.add_argument_group('advanced')
    adv.add_argument('--trim-front', type=int, default=0, metavar='N',
                     help='remove N bases from the start of each sequence before analysis')
    adv.add_argument('--trim-end', type=int, default=0, metavar='N',
                     help='remove N bases from the end of each sequence before analysis')

    p.add_argument('--visual', action='store_true',
                   help='also print a monospace alignment of each stretch')
    p.add_argument('--matrix', action='store_true',
                   help='also print pairwise Hamming distances between motifs found')
    p.add_argument('--expand', action='store_true',
                   help='treat inputs as nomenclature and print the DNA sequence')
    p.add_argument('--color', choices=('auto', 'always', 'never'), default='auto',
                   help='ANSI colour in --visual / --matrix output (default: auto)')
    args = p.parse_args()

    inputs = gather_inputs(args)
    if not inputs:
        p.error('no input given (pass sequences, --input FILE, or pipe via stdin)')

    if args.expand:
        for line in inputs:
            print(convert_nomenclature_to_sequence(line.strip().upper()))
        return

    palette = Palette(args.color == 'always'
                      or (args.color == 'auto' and sys.stdout.isatty()))
    structured = args.visual or args.matrix
    all_canonicals = set()

    for idx, line in enumerate(inputs, 1):
        seq, had_invalid = prepare_sequence(line)
        if args.trim_front or args.trim_end:
            end = len(seq) - args.trim_end if args.trim_end else len(seq)
            seq = seq[args.trim_front:max(args.trim_front, end)]
        if args.orientation == 'reverse-complement':
            seq = reverse_complement(seq)
        stretches = call_stretches(seq, args.min_repeats)
        for s in stretches:
            all_canonicals.add(s['canonical'])

        color_index = {c: i for i, c in
                       enumerate(sorted({s['canonical'] for s in stretches}))}
        nomenclature, segments = build_nomenclature(
            seq, stretches, args.hide_n, args.full_seq_gaps, args.strip_ends_n)

        if structured:
            if idx > 1:
                print()
            print(f"# Result {idx}")
            if had_invalid:
                print("# (invalid characters were ignored)")
            print(f"seq  : {seq}")
            print(f"nom  : {nomenclature or '(no repeats found)'}")
            if args.visual and segments:
                print(build_visual(seq, segments, color_index, palette))
        else:
            print(nomenclature or '(no repeats found)')

    if args.matrix:
        matrix = build_matrices(all_canonicals, palette)
        if matrix:
            print()
            print(matrix)


if __name__ == '__main__':
    main()