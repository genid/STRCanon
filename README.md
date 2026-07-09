# STRCanon

[![Open in browser](https://img.shields.io/badge/Open%20web%20app-4a7ff5?style=for-the-badge&logo=googlechrome&logoColor=white)](https://genid.github.io/STRCanon/)

A simple, unbiased canonical motif based system for short tandem repeat (STR) nomenclature.

Short tandem repeats can be described starting at any rotation of their
repeat unit, and on either DNA strand, so the same repeat is often written
several different ways (e.g. `AAG`, `AGA`, `GAA`, or their reverse
complements `CTT`, `TTC`, `TCT`). STRCanon removes that ambiguity by mapping
every rotation/strand variant of a motif to a single canonical form, then
renders sequences in a compact bracket nomenclature built on top of it.

## Canonicalization rule

For every k-mer (k = 1..6), build its equivalence class from all cyclic
rotations of the motif and all cyclic rotations of its reverse complement.
The alphabetically-first sequence in that class is chosen as the canonical
representative. Motifs that are themselves a repeat of a shorter unit (e.g.
`ATAT`, a repeat of `AT`) are excluded, since they aren't valid independent
motifs at that length.

## Nomenclature

```
[AAG]        canonical motif, no rotation
[AAG>1]      canonical rotated 1 position left (AAG -> AGA)
[~AAG]       reverse complement of canonical (AAG -> CTT)
[~AAG>2]     RC of canonical rotated 2 left  (AAG -> GAA -> TTC)
]3           3 complete repeats
]3.1         3 complete repeats + 1 partial base
[N]5         5-base gap between stretches
[N]-2        2-base overlap between stretches
(ACGT)       literal gap sequence (with --full-seq-gaps)
```

A stretch may claim a partial repeat whose bases the *next* stretch also
covers — `[AAG]3.1` and `[AG]4` in the example below both describe the `A`
at offset 9. `[N]-2` records exactly that: the following stretch starts 2
bases before the previous one ended. When expanding, terms are therefore
written through a cursor rather than concatenated: `[N]-2` rewinds the
cursor 2 bases and the next stretch overwrites them.

```
AAGAAGAAGAGAGAG   ->  [AAG]3.1[N]-3[AG]4
AAGAAGAAGA             [AAG]3.1
       AGAGAGAG                 [AG]4   (starts 3 bases early -> [N]-3)
```

Counts use plain digits, so any output line can be fed straight back in with
`--expand`. The reconstruction is exact when no information was suppressed on
the way out — that is, with `--no-strip-ends-n` (keep the flanking sequence)
and the default `--full-seq-gaps` (write gap bases literally). Under
`--no-full-seq-gaps` a gap becomes `[N]n`, which expands to `n` Ns: the length
is preserved but the gap bases are not recoverable. Under the default
`--strip-ends-n` the flanks are dropped, so `--expand` returns the STR region
only.

## Files

- **`strcanon.py`** — command-line tool. Detects STR stretches in DNA
  sequences and renders them in canonical bracket nomenclature. The
  canonical motif for every k-mer is derived on the fly using the same
  lexicographic-minimum rule as `generate_canonicals.py`, so no external
  lookup table is required at run time.
- **`generate_canonicals.py`** — generates the canonical motif lookup
  table for all k-mers of length 1-6 and writes it to
  `str_canonical_motifs.tsv` (columns: `motif`, `canonical`,
  `is_reverse_complement`, `k`).
- **`str_canonical_motifs.tsv`** — precomputed lookup table produced by
  `generate_canonicals.py`, provided for reference and for use in other
  tools/languages.
- **`test_strcanon.py`** — test suite (`python -m unittest test_strcanon`).
  Covers canonicalization against the lookup table, rendering of common
  forensic loci (D21S11, vWA, FGA, TH01, …), strand symmetry, and
  `expand(render(seq)) == seq` over randomized repeat-dense sequences.
- **`index.html`** — a self-contained, in-browser version of the tool
  ([live at genid.github.io/STRCanon](https://genid.github.io/STRCanon/)):
  paste in sequences or nomenclature strings and get canonical nomenclature,
  visual alignments, and pairwise motif distance matrices, with no
  installation required. Open the file directly in a browser.

## Command-line usage

```bash
strcanon.py AAGAAGAAGAGAGAG
strcanon.py TGCAAGAAGAAG --no-strip-ends-n          # keep the leading (TGC)
strcanon.py --input seqs.txt --visual --matrix
echo TTCTTCTTCTTCTT | strcanon.py                   # -> [~AAG>2]4.2
strcanon.py --expand '(TGC)[AAG]4[AG]3.1'           # nomenclature -> sequence
```

Options:

| Flag | Description |
| --- | --- |
| `sequences` | DNA sequence(s) or nomenclature string(s) to analyse |
| `-i, --input FILE` | read one sequence/nomenclature per line from FILE |
| `--min-repeats N` | minimum repeats to call a stretch (mono-nucleotides always require 5; default: 3) |
| `--hide-n / --no-hide-n` | suppress all `[N]` gap/overlap markers |
| `--full-seq-gaps / --no-full-seq-gaps` | show literal bases in `(parentheses)` instead of `[N]n` (default: on) |
| `--strip-ends-n / --no-strip-ends-n` | suppress only the leading/trailing gap markers (default: on) |
| `--trim-front N` | remove N bases from the start of each sequence before analysis |
| `--trim-end N` | remove N bases from the end of each sequence before analysis |
| `--visual` | also print a monospace alignment of each stretch |
| `--matrix` | also print pairwise Hamming distances between motifs found |
| `--expand` | treat inputs as nomenclature and print the DNA sequence |
| `--color {auto,always,never}` | ANSI colour in `--visual` / `--matrix` output (default: auto) |

## Requirements

Python 3.9+ (standard library only, no dependencies).

## Tests

```bash
python -m unittest test_strcanon
```

## Generating the lookup table

```bash
python generate_canonicals.py
```

This regenerates `str_canonical_motifs.tsv` from scratch.
