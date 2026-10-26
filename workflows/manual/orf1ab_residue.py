#!/usr/bin/env python3
"""
Work out the true amino-acid consequence of an ORF1ab variant, by building the
CDS from the GFF join and translating it. Independent of snpEff.

It prints the answer two ways: with the slippage nucleotide kept (correct, as
GenBank annotates it) and with it collapsed (what snpEff appears to have done),
so you can see exactly where the two disagree.

Usage:
    python3 orf1ab_residue.py <ref.fasta> <ref.gff> <position> <ALT_base>
Example:
    python3 orf1ab_residue.py reference/MZ558051.1.fasta reference/MZ558051.1.gff 13977 T
"""

import sys

CODONS = {}
_bases = "TCAG"
_aa = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
_i = 0
for b1 in _bases:
    for b2 in _bases:
        for b3 in _bases:
            CODONS[b1 + b2 + b3] = _aa[_i]
            _i += 1

THREE = {"A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
         "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
         "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
         "Y": "Tyr", "V": "Val", "*": "Ter"}


def read_fasta(path):
    """Return the single sequence in a FASTA file, uppercased, as one string."""
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def orf1ab_segments(gff_path):
    """The CDS segments of the ORF1ab polyprotein (the longest ORF1ab protein).

    A GFF may hold both ORF1ab and the shorter ORF1a under gene=ORF1ab, so we
    group by protein_id and keep whichever is longest.
    """
    by_protein = {}
    for line in open(gff_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "CDS" or "gene=ORF1ab" not in f[8]:
            continue
        pid = None
        for kv in f[8].split(";"):
            if kv.startswith("protein_id="):
                pid = kv.split("=", 1)[1]
        by_protein.setdefault(pid, []).append((int(f[3]), int(f[4])))
    best = max(by_protein.items(), key=lambda kv: sum(e - s + 1 for s, e in kv[1]))
    pid, segs = best
    return pid, sorted(segs)


def build_cds(seq, segments, keep_duplicate):
    """Concatenate the segments. Returns the CDS string and a map from each
    CDS index (0-based) back to its genomic position.

    keep_duplicate=True  -> the slippage base appears twice (correct)
    keep_duplicate=False -> the second segment starts one base later
    """
    cds, gmap = [], []
    previous_end = None
    for start, end in segments:
        s = start
        if not keep_duplicate and previous_end is not None and start <= previous_end:
            s = previous_end + 1
        for g in range(s, end + 1):
            cds.append(seq[g - 1])
            gmap.append(g)
        previous_end = end
    return "".join(cds), gmap


def report(label, cds, gmap, position, alt):
    if position not in gmap:
        print(f"{label}: position {position} is not in this CDS model")
        return
    idx = gmap.index(position)                 # 0-based CDS index
    codon_no = idx // 3 + 1                    # 1-based residue number
    frame = idx % 3                            # 0,1,2 within the codon
    start = idx - frame
    ref_codon = cds[start:start + 3]
    mut_codon = ref_codon[:frame] + alt + ref_codon[frame + 1:]
    ref_aa = CODONS.get(ref_codon, "?")
    mut_aa = CODONS.get(mut_codon, "?")
    genomic = [gmap[start + i] for i in range(3)]
    kind = "synonymous" if ref_aa == mut_aa else ("nonsense" if mut_aa == "*" else "missense")
    print(f"{label}")
    print(f"    CDS length           {len(cds)} nt = {len(cds)//3} codons ({len(cds)//3 - 1} aa + stop)")
    print(f"    CDS position         c.{idx+1}{cds[idx]}>{alt}   (base {frame+1} of the codon)")
    print(f"    codon {codon_no:<5}          genomic {genomic[0]}-{genomic[2]}  {ref_codon} -> {mut_codon}")
    print(f"    consequence          p.{THREE.get(ref_aa,'?')}{codon_no}{THREE.get(mut_aa,'?')}  ({kind})")


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    fasta, gff, position, alt = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4].upper()
    seq = read_fasta(fasta)
    pid, segs = orf1ab_segments(gff)

    print(f"ORF1ab protein {pid}, segments " + ", ".join(f"{s}..{e}" for s, e in segs))
    print(f"reference base at genomic {position}: {seq[position-1]}  (variant -> {alt})")
    print()
    cds, gmap = build_cds(seq, segs, keep_duplicate=True)
    report("CORRECT - slippage base kept (as GenBank annotates it):", cds, gmap, position, alt)
    print()
    cds2, gmap2 = build_cds(seq, segs, keep_duplicate=False)
    report("WRONG - slippage base collapsed (what snpEff appears to have done):",
           cds2, gmap2, position, alt)


if __name__ == "__main__":
    main()
