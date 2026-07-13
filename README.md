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

## Analysis strand

Repeats are detected 5'→3': a stretch is matched left to right, and a partial
repeat is whatever is left over on the run's 3' end. Reverse-complementing a
sequence moves that leftover to the 5' end, so calling stretches on whichever
strand happened to be sequenced would place partial repeats differently. TH01
allele 9.3 is the textbook case — its 3-base `ATG` interruption gets absorbed
into the first stretch on one strand and the second on the other.

So stretches are always called on **one strand of the duplex**, picked from the
molecule rather than from the input: the lexicographically smaller of the
sequence and its reverse complement. Both strands of a locus canonicalize to
the same string, so both yield the same blocks. The result is then mirrored
onto whichever strand is being displayed.

`--orientation` therefore changes only the order the motifs are shown in and
how each is written (`~` and `>r`) — never the canonical motif, the repeat
count, or the partial length:

```
TH01 allele 9.3, plus strand   [AATG]6.1[AATG>2]3.2
             ... minus strand  [~AATG]3.2[~AATG>1]6.1
```

Both name the same two blocks — `AATG` × 6 + 1, and `AATG` × 3 + 2 — in
mirror-image order.

Because the analysis strand is chosen by comparing the whole sequence,
extending the flanks far enough to flip that comparison can move a partial
repeat from one end of a run to the other. Trim consistently (`--trim-front` /
`--trim-end`) when comparing alleles across samples.

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
{TH01 CE=9.3 region=42bp len=125bp}
             CE annotation prefix (see below); not part of the repeat
             grammar, so --expand ignores it
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

## CE alleles and forensic markers

A capillary electrophoresis (CE) allele is a *length* call: the fragment is
sized, and its size is read off as a repeat count. Sequence nomenclature carries
far more information, but it does not by itself say which CE allele a sequence
would have been called as — and CE allele numbers are what the rest of forensic
practice is written in (databases, match reports, published frequencies). So
STRCanon reports both:

```
$ strcanon.py TGCAGGTCACAGGGAACACAGACTCCATGGTGAATGAATG...GGGAAATAAGGGAGGAACAGG
{TH01 CE=9.3 region=42bp len=125bp} [AATG>2]6.2[AATG>1]4
```

The pair is the point. Two sequences with the **same CE allele but different
brackets** are iso-alleles: one length, two molecules, which CE alone cannot
tell apart. Real example — both of these are TH01 allele 9:

```
{TH01 CE=9 region=39bp len=230bp} [AATG>2]9.3      MH085122.2
{TH01 CE=9 region=39bp len=230bp} [AATG>2]8.3      MH085123.2   (TGAG first repeat)
```

The sequence length is reported whether or not a marker is recognised. The
braces sit outside the repeat grammar, so an annotated line still round-trips
through `--expand`; `--no-ce` turns the prefix off entirely.

### How the CE allele is inferred

Each marker in the panel supplies two flanking **anchor** sequences, the repeat
**period** (4 = tetranucleotide, 5 = penta, …), and an **offset** — the number
of non-repeat bases between the anchors, i.e. the region length of a
hypothetical allele 0. The region between the anchors is then

```
region_len = offset + period * allele
```

so the allele is `divmod(region_len - offset, period)`: the quotient is the
repeat count and the remainder is the `.n` of a microvariant. TH01 9.3 is 9
whole `AATG` repeats plus a 3-base `ATG`, and lands 42 bases between its
anchors: `3 + 4×9 + 3`.

The anchors are what make this kit-independent. An amplicon's absolute size
depends on where a kit's primers sit, but the offset is calibrated against the
*same* anchor pair used to measure the region, so the primer-dependent constant
cancels out: any sequence containing both anchors gives the same CE allele, no
matter how much flanking sequence came with it.

Anchors tolerate one mismatch by default (`--flank-mismatches`), because
flanking SNPs are common — STRSeq MH085118.2 is a TH01 allele 7 carrying
`AACAGAGACT` where the reference reads `AACACAGACT`, and an exact-only anchor
misses the locus entirely.

### When there is no call

The tool declines to guess in two situations, and says which:

- **No flanking sequence.** A sequence trimmed down to the repeat array itself
  has nothing to anchor to, so its length is not recoverable and no CE allele is
  reported. It is still named — just not sized.
- **Ambiguous anchors.** If an anchor matches in more than one place, the STR
  region has no single length. STRSeq MN983127.1 is a D13S317 duplication allele
  (28.2) in which the whole repeat block *and its 3′ flank* occur twice, so the
  3′ anchor matches once *inside* the array. Taking the first match would report
  a confident allele 10 for a sequence that is really a 28.2, so STRCanon reports
  `CE=ambiguous` instead.

### The panel

The bundled panel is **47 loci — 25 autosomal, 5 X-STR and 17 Y-STR** — covering
the CODIS core, Penta D/E, and the X and Y markers of the common kits. Its
anchors and offsets were fitted to [NIST STRSeq](https://strseq.nist.gov/)
(BioProject PRJNA380127), the public catalogue of sequenced STR alleles.

Each anchor pair is pulled in as close to the repeat as it can go while staying
unambiguous, and is kept only if

    offset = region_len - (period × repeats + remainder)

comes out **identical for every record it covers**. A pair that brackets the
variable region and nothing else gives one offset; any other pair gives several,
and is rejected. Nothing is taken on trust — a locus is admitted only if its own
records agree.

Run blind over all 2578 STRSeq records — the panel is given the sequence and
nothing else, and must pick the locus itself. Of those, 1836 belong to a locus
the panel carries:

| | |
| --- | --- |
| called | **1555** — every one at the **right locus** with the **published allele** |
| declined | 280 unanchorable (too little flank), 1 ambiguous (the duplication) |
| mis-assigned | **0** |
| wrong allele | **0** |

Nothing is called wrong; what cannot be determined is declined. That matters
more once X and Y markers are in the panel, because several of their anchors are
themselves repeat-like (DYS391's is `TATCTGTCTGTCTG`) and could in principle
match inside another locus's array — the blind run is what shows they do not.

The unanchorable records are a property of the data, not the method — STRSeq
records are trimmed to the ISFG reported range, and many retain under 12 bp of
flank. The panel's anchors therefore sit as close to the repeat as they can
while staying unambiguous; anchors placed further out (as read-based tools use,
since they see whole amplicons) fall outside that range and are simply absent.

Some loci are deliberately **not** in the panel. SE33, DXS10135, DYS458, D4S2408
and Y-GATA-H4 derive a consistent offset but their anchors reach only a small
fraction of their own records, so they would almost never fire; DYS385a/b is a
two-copy locus, for which a single length call is not well defined.

Panels are swappable. `--markers FILE` takes the same 7-column layout, which is
also STRait Razor's locus-config format, so a lab's existing kit config —
ForenSeq, PowerSeq, SE33, anything not in the bundled panel — can be passed
in directly:

```bash
strcanon.py --list-markers                      # show the bundled panel
strcanon.py --markers ForenSeq.config seqs.txt  # use a lab panel instead
strcanon.py --marker TH01 --input seqs.txt      # force one marker
```

**Validate any panel against your own kit before casework.**

## Files

- **`strcanon.py`** — command-line tool. Detects STR stretches in DNA
  sequences and renders them in canonical bracket nomenclature. The
  canonical motif for every k-mer is derived on the fly using the same
  lexicographic-minimum rule as `generate_canonicals.py`, so no external
  lookup table is required at run time.
- **`str_markers.py`** — the forensic marker panel and CE allele inference:
  anchor matching, the `offset + period × allele` arithmetic, and the panel
  itself. Run it directly to print the panel; `--export` regenerates
  `str_markers.tsv`.
- **`str_markers.tsv`** — the panel exported for use by other tools/languages
  (columns: `marker`, `type`, `flank5`, `flank3`, `motif`, `period`, `offset`).
- **`generate_canonicals.py`** — generates the canonical motif lookup
  table for all k-mers of length 1-6 and writes it to
  `str_canonical_motifs.tsv` (columns: `motif`, `canonical`,
  `is_reverse_complement`, `k`).
- **`str_canonical_motifs.tsv`** — precomputed lookup table produced by
  `generate_canonicals.py`, provided for reference and for use in other
  tools/languages.
- **`test_strcanon.py`** — test suite (`python -m unittest test_strcanon`).
  Covers canonicalization against the lookup table, rendering of common
  forensic loci (D21S11, vWA, FGA, TH01, …), strand symmetry of the calls, and
  `expand(render(seq)) == seq` on both strands over randomized repeat-dense
  sequences.
- **`test_markers.py`** — test suite for the CE allele feature
  (`python -m unittest test_markers`). The fixtures are real STRSeq records
  carrying their published length-based allele, so the panel is regression-locked
  against ground truth: the TH01 9.3 microvariant, a flanking-SNP allele that is
  only found because anchors tolerate a mismatch, a record trimmed past its
  anchors (no call), the D13S317 duplication (ambiguous), and an iso-allele pair.
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
strcanon.py --input alleles.txt --no-ce             # sequence nomenclature only
strcanon.py --list-markers                          # the forensic marker panel
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
| `--orientation {as-given,reverse-complement}` | strand to report on; affects only display order and how each motif is written, never how it is called (default: as-given) |
| `--trim-front N` | remove N bases from the start of each sequence before analysis |
| `--trim-end N` | remove N bases from the end of each sequence before analysis |
| `--visual` | also print a monospace alignment of each stretch |
| `--matrix` | also print pairwise Hamming distances between motifs found |
| `--expand` | treat inputs as nomenclature and print the DNA sequence |
| `--color {auto,always,never}` | ANSI colour in `--visual` / `--matrix` output (default: auto) |
| `--ce / --no-ce` | prefix each line with the sequence length and, when a panel marker is recognised, its CE allele (default: on) |
| `--markers FILE` | marker panel to use instead of the bundled one (also reads STRait Razor `.config` files) |
| `--marker NAME` | force a single marker instead of searching the whole panel |
| `--flank-mismatches N` | mismatches tolerated per anchor, for flanking SNPs (default: 1) |
| `--list-markers` | print the marker panel and exit |

## Requirements

Python 3.9+ (standard library only, no dependencies).

## Tests

```bash
python -m unittest discover -p "test_*.py"
```

## Generating the lookup table

```bash
python generate_canonicals.py
```

This regenerates `str_canonical_motifs.tsv` from scratch.
