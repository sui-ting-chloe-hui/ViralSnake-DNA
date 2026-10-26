#!/usr/bin/env python3
"""
Independently re-derive the gene and codon number for a list of genomic
positions, straight from the reference GFF3. Nothing here depends on snpEff,
so it is a genuine cross-check rather than a restatement.

Usage:
    python3 verify_annotation.py <reference.gff3>
"""

import sys
from collections import defaultdict

# The positions our report annotates, and what it claims for each.
# Format: genomic position -> (gene, codon number claimed in the report)
REPORT_CLAIMS = {
    5168:  ("ORF1ab", 1640),
    11067: ("ORF1ab", 3606),
    13977: ("ORF1ab", 4576),
    21990: ("S",      148),
    21992: ("S",      149),
    28852: ("N",      199),
    29653: ("ORF10",  38),
}

# Wuhan-Hu-1 (NC_045512.2) CDS starts, for showing the coordinate offset.
WUHAN_STARTS = {"ORF1ab": 266, "S": 21563, "N": 28274, "ORF10": 29558}


def read_cds(gff_path):
    """Collect CDS intervals per gene, in genomic order.

    A gene can have several CDS lines. ORF1ab has two because of the -1
    ribosomal frameshift: GenBank writes it as join(a..b, b..c), where the
    nucleotide at b is deliberately listed twice. Keeping both segments is
    what makes downstream codon numbering correct.
    """
    cds = defaultdict(list)
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            start, end, attrs = int(f[3]), int(f[4]), f[8]
            gene = None
            for kv in attrs.split(";"):
                if kv.startswith("gene="):
                    gene = kv.split("=", 1)[1]
                    break
            if gene:
                cds[gene].append((start, end))
    for gene in cds:
        cds[gene].sort()
    return cds


def codon_number(position, segments):
    """Return the 1-based codon (amino-acid) number for a genomic position.

    We walk the CDS segments in order, accumulating length, until we find the
    segment containing the position. Because the slippage nucleotide appears
    in two segments, positions after the slip site are correctly shifted by
    one relative to a naive end-to-end calculation.
    """
    consumed = 0
    for start, end in segments:
        if start <= position <= end:
            cds_pos = consumed + (position - start + 1)
            return -(-cds_pos // 3), cds_pos   # ceiling division
        consumed += end - start + 1
    return None, None


def main():
    if len(sys.argv) != 2:
        sys.exit("usage: python3 verify_annotation.py <reference.gff3>")
    cds = read_cds(sys.argv[1])

    print("Gene boundaries as annotated in this GFF")
    print("-" * 78)
    for gene in sorted(cds):
        segs = ", ".join(f"{s}..{e}" for s, e in cds[gene])
        note = ""
        if gene in WUHAN_STARTS:
            note = f"   (Wuhan-Hu-1 start {WUHAN_STARTS[gene]}, offset {cds[gene][0][0] - WUHAN_STARTS[gene]:+d})"
        print(f"  {gene:<10} {segs}{note}")

    print()
    print("Re-derived codon numbers vs. the report")
    print("-" * 78)
    print(f"{'position':>9}  {'gene':<8} {'CDS pos':>8} {'derived':>8} {'report':>7}  verdict")
    failures = 0
    for pos in sorted(REPORT_CLAIMS):
        claim_gene, claim_codon = REPORT_CLAIMS[pos]
        hit_gene, derived, cds_pos = None, None, None
        for gene, segs in cds.items():
            c, p = codon_number(pos, segs)
            if c is not None:
                hit_gene, derived, cds_pos = gene, c, p
                break
        ok = (hit_gene == claim_gene and derived == claim_codon)
        if not ok:
            failures += 1
        print(f"{pos:>9}  {str(hit_gene):<8} {str(cds_pos):>8} {str(derived):>8} {claim_codon:>7}  "
              f"{'match' if ok else 'MISMATCH'}")

    print()
    print(f"{len(REPORT_CLAIMS) - failures} of {len(REPORT_CLAIMS)} match; {failures} mismatch.")


if __name__ == "__main__":
    main()
