#!/usr/bin/env python3
"""Generate canonical short tandem repeat (STR) motifs for k-mers of length 1-6.

Canonicalization rule: within each equivalence class of a motif (all cyclic
rotations, plus all cyclic rotations of its reverse complement), the
alphabetically-first sequence is chosen as the canonical representative.

Output columns:
  motif      - the k-mer (rotation or reverse complement variant)
  canonical  - the canonical base motif it maps to
  is_rc      - 1 if this motif is a reverse complement of the canonical, 0 otherwise
  k          - length of the motif
"""

from itertools import product

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}


def reverse_complement(seq: str) -> str:
    return ''.join(COMPLEMENT[b] for b in reversed(seq))


def get_rotations(seq: str) -> list[str]:
    n = len(seq)
    return [seq[i:] + seq[:i] for i in range(n)]


def has_sub_period(seq: str) -> bool:
    """Return True if seq is a tandem repeat of a shorter unit."""
    n = len(seq)
    for period in range(1, n):
        if n % period == 0 and seq == seq[:period] * (n // period):
            return True
    return False


def generate_lookup(k_max: int = 6) -> tuple[dict, list[str]]:
    """Build the motif -> (canonical, is_rc) lookup table for k = 1..k_max."""
    lookup: dict = {}
    canonicals: list[str] = []

    for k in range(1, k_max + 1):
        kmers = sorted(''.join(p) for p in product('ACGT', repeat=k))
        processed: set[str] = set()

        for kmer in kmers:
            if kmer in processed:
                continue

            if has_sub_period(kmer):
                # Degenerate motif (a repeat of a shorter unit) -- not a
                # valid independent STR motif at this length.
                for rot in get_rotations(kmer):
                    lookup.setdefault(rot, None)
                    lookup.setdefault(reverse_complement(rot), None)
                    processed.add(rot)
                    processed.add(reverse_complement(rot))
                continue

            fwd = set(get_rotations(kmer))
            rc = set(get_rotations(reverse_complement(kmer)))
            group = fwd | rc

            canonical = sorted(group)[0]
            canonicals.append(canonical)

            for motif in group:
                is_rc = motif not in fwd
                lookup[motif] = (canonical, is_rc)

            processed.update(group)

    return lookup, canonicals


def write_tsv(lookup: dict, output_file: str = 'str_canonical_motifs.tsv') -> None:
    entries = [
        (motif, canonical, is_rc)
        for motif, val in lookup.items()
        if val is not None
        for canonical, is_rc in (val,)
    ]
    entries.sort(key=lambda x: (len(x[0]), x[0]))

    with open(output_file, 'w') as f:
        f.write('motif\tcanonical\tis_reverse_complement\tk\n')
        for motif, canonical, is_rc in entries:
            f.write(f'{motif}\t{canonical}\t{1 if is_rc else 0}\t{len(motif)}\n')

    print(f"Lookup table written to: {output_file}")
    print(f"Total entries: {len(entries)}")


def main() -> None:
    lookup, canonicals = generate_lookup(k_max=6)
    write_tsv(lookup)

    print()
    print("Example check (AGAT group):")
    for motif in ['AGAT', 'GATA', 'ATAG', 'TAGA', 'ATCT', 'TCTA', 'CTAT', 'TATC']:
        val = lookup.get(motif)
        if val:
            canonical, is_rc = val
            print(f"  {motif} -> canonical={canonical}, is_rc={int(is_rc)}")


if __name__ == '__main__':
    main()