#!/usr/bin/env python3
"""
Build the SARS-CoV-2 variant report from ViralSnake-DNA pipeline output.

Reads, for each sample, from <proj_dir>/results/02.variants/:
  <sample>.combined.annotate.vcf           GATK high-confidence SNPs+indels, snpEff-annotated
  <sample>.minor_variants.reheader.vcf.gz  LoFreq minor variants (falls back to .vcf / non-reheader)

Writes:
  <out>.docx   the report, with the same sections and tables as the manual version
  <out>.tsv    every variant in one tidy table, for checking and for the shared drive

Usage:
  python variant_report.py --proj-dir projects/SARS2-BavPat1 \
                           --samples p10CCT1 p10CCT2 p10CCT3 \
                           --out sars2.report
  # --samples omitted: every sample with a combined.annotate.vcf is picked up
"""

import argparse
import csv
import glob
import gzip
import os
import re
import sys
from collections import Counter

# ---------------------------------------------------------------- VCF parsing

# snpEff ANN subfields, in order, as documented in the ANN spec.
ANN_FIELDS = [
    "allele", "annotation", "impact", "gene_name", "gene_id", "feature_type",
    "feature_id", "transcript_biotype", "rank", "hgvs_c", "hgvs_p",
    "cdna_pos", "cds_pos", "aa_pos", "distance", "errors",
]

# snpEff annotation terms we treat as regulatory rather than protein-altering.
REGULATORY_TERMS = {
    "upstream_gene_variant", "downstream_gene_variant",
    "intergenic_region", "intragenic_variant",
    "5_prime_UTR_variant", "3_prime_UTR_variant",
    "5_prime_UTR_premature_start_codon_gain_variant",
}

THREE_LETTER = {
    "Ala": "Ala", "Arg": "Arg", "Asn": "Asn", "Asp": "Asp", "Cys": "Cys",
    "Gln": "Gln", "Glu": "Glu", "Gly": "Gly", "His": "His", "Ile": "Ile",
    "Leu": "Leu", "Lys": "Lys", "Met": "Met", "Phe": "Phe", "Pro": "Pro",
    "Ser": "Ser", "Thr": "Thr", "Trp": "Trp", "Tyr": "Tyr", "Val": "Val",
}


def open_maybe_gzip(path):
    """Open a plain or gzipped text file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_info(info_string):
    """Turn an INFO column into a dict. Flags with no '=' get the value True."""
    info = {}
    for field in info_string.split(";"):
        if not field:
            continue
        if "=" in field:
            key, value = field.split("=", 1)
            info[key] = value
        else:
            info[field] = True
    return info


def best_annotation(info):
    """
    Return the first (highest-ranked) snpEff ANN entry as a dict.
    snpEff orders ANN by predicted severity, so entry 0 is the one to report.
    """
    raw = info.get("ANN")
    if not raw:
        return {}
    first = raw.split(",")[0].split("|")
    return dict(zip(ANN_FIELDS, first))


def tidy_protein_change(hgvs_p):
    """'p.Asp420Tyr' -> 'Asp420Tyr'. Returns '' when there is no protein change."""
    if not hgvs_p:
        return ""
    return re.sub(r"^p\.", "", hgvs_p)


def classify(ref, alt, annotation, impact):
    """
    Sort a variant into one of the report's three categories.
    Indels are decided by length; SNPs by the snpEff annotation term.
    """
    if len(ref) != len(alt):
        return "indel"
    if annotation in REGULATORY_TERMS:
        return "regulatory"
    if impact in ("HIGH", "MODERATE"):
        return "protein_altering"
    # synonymous_variant and anything else coding but silent
    return "other"


def indel_type(ref, alt):
    """Human-readable indel label, matching the wording in the manual report."""
    if len(alt) > len(ref):
        inserted = len(alt) - len(ref)
        return "Large Insertion" if inserted >= 10 else "Insertion"
    return "Deletion"


def shorten(seq, keep=6):
    """Trim a long allele for display: ACCTCGGTTA... -> ACCTCG..."""
    if len(seq) <= keep + 3:
        return seq
    return seq[:keep] + "..."


def read_annotated_vcf(path, sample, pass_only=True):
    """Read one <sample>.combined.annotate.vcf into a list of dicts."""
    variants = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info_string = col[:8]
            if pass_only and filt not in ("PASS", ".", ""):
                continue
            # A multi-allelic record would need splitting; flag it rather than guess.
            if "," in alt:
                print(f"  ! {sample} {chrom}:{pos} is multi-allelic ({alt}) - "
                      f"reporting the first ALT only", file=sys.stderr)
                alt = alt.split(",")[0]
            info = parse_info(info_string)
            ann = best_annotation(info)
            annotation = ann.get("annotation", "")
            impact = ann.get("impact", "")
            variants.append({
                "sample": sample,
                "source": "GATK",
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filt,
                "gene": ann.get("gene_name", ""),
                "annotation": annotation,
                "impact": impact,
                "protein_change": tidy_protein_change(ann.get("hgvs_p", "")),
                "category": classify(ref, alt, annotation, impact),
                "depth": info.get("DP", ""),
                "af": "",
            })
    return variants


def read_minor_vcf(path, sample, min_af, pass_only=True):
    """Read one LoFreq minor-variant VCF. AF and DP come from the INFO column."""
    variants = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            col = line.rstrip("\n").split("\t")
            if len(col) < 8:
                continue
            chrom, pos, _id, ref, alt, qual, filt, info_string = col[:8]
            if pass_only and filt not in ("PASS", ".", ""):
                continue
            info = parse_info(info_string)
            try:
                af = float(info.get("AF", "nan"))
            except ValueError:
                af = float("nan")
            if af != af or af < min_af:          # af != af is the NaN test
                continue
            ann = best_annotation(info)
            variants.append({
                "sample": sample,
                "source": "LoFreq",
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "qual": qual,
                "filter": filt,
                "gene": ann.get("gene_name", ""),
                "annotation": ann.get("annotation", ""),
                "impact": ann.get("impact", ""),
                "protein_change": tidy_protein_change(ann.get("hgvs_p", "")),
                "category": "minor",
                "depth": info.get("DP", ""),
                "af": f"{af:.4f}",
            })
    return variants


def find_minor_vcf(variants_dir, sample):
    """
    The minor-variant file has had several names across pipeline versions.
    Try the reheadered one first, then the plain one, then the annotated LoFreq VCF.
    """
    candidates = [
        f"{sample}.minor_variants.reheader.vcf.gz",
        f"{sample}.minor_variants.reheader.vcf",
        f"{sample}.minor_variants.vcf.gz",
        f"{sample}.minor_variants.vcf",
        f"{sample}.lofreq.filtered.ann.vcf",
    ]
    for name in candidates:
        path = os.path.join(variants_dir, name)
        if os.path.exists(path):
            return path
    return None


# ------------------------------------------------------------------- document

def add_table(document, headers, rows, bold_predicate=None):
    """
    Add a table with a header row. bold_predicate(row) -> True bolds that row,
    which is how the fixed (>= fixed_af) minor variants are marked.
    """
    table = document.add_table(rows=1, cols=len(headers))
    table.style = "Table Grid"
    for cell, text in zip(table.rows[0].cells, headers):
        cell.text = ""
        run = cell.paragraphs[0].add_run(str(text))
        run.bold = True
    for row in rows:
        cells = table.add_row().cells
        emphasise = bool(bold_predicate and bold_predicate(row))
        for cell, value in zip(cells, row):
            cell.text = ""
            run = cell.paragraphs[0].add_run("" if value is None else str(value))
            run.bold = emphasise
    document.add_paragraph()
    return table


def group_samples(rows, key_columns):
    """
    Collapse rows that are identical except for the sample, so a variant shared
    by two samples becomes one row reading 'p10CCT1, p10CCT2' as in the manual report.
    """
    grouped = {}
    order = []
    for row in rows:
        key = tuple(row[column] for column in key_columns)
        if key not in grouped:
            grouped[key] = []
            order.append(key)
        grouped[key].append(row["sample"])
    return [(", ".join(sorted(set(samples))), key) for key in order
            for samples in [grouped[key]]]


def build_document(samples, all_variants, args):
    from docx import Document
    from docx.shared import Pt

    document = Document()
    if args.title:
        document.add_heading(args.title, level=0)

    gatk = [v for v in all_variants if v["source"] == "GATK"]
    minor = [v for v in all_variants if v["source"] == "LoFreq"]

    protein = [v for v in gatk if v["category"] == "protein_altering"]
    regulatory = [v for v in gatk if v["category"] == "regulatory"]
    indels = [v for v in gatk if v["category"] == "indel"]
    other = [v for v in gatk if v["category"] == "other"]

    # --- 1. Key summary ---------------------------------------------------
    document.add_heading("1. Key Summary - SNP and INDEL analysis", level=1)
    sample_list = ", ".join(samples)
    document.add_paragraph(
        f"This report summarizes the genetic variants identified in samples "
        f"{sample_list}, called against {args.reference}. This analysis identified "
        f"high-confidence genetic variants using GATK's strict filtering criteria. "
        f"The variants are categorised by their functional impact:"
    )

    def gene_note(variants):
        genes = sorted({v["gene"] for v in variants if v["gene"]})
        return f" ({', '.join(genes)} gene{'s' if len(genes) > 1 else ''})" if genes else ""

    document.add_paragraph(
        f"{len(protein)} Protein-Altering SNP{'s' if len(protein) != 1 else ''}:"
        f" coding changes{gene_note(protein)}"
    )
    document.add_paragraph(
        f"{len(regulatory)} Regulatory SNP{'s' if len(regulatory) != 1 else ''}:"
        f" upstream/downstream variants affecting regulatory regions"
    )
    document.add_paragraph(
        f"{len(indels)} Insertion/Deletion{'s' if len(indels) != 1 else ''} (Indels):"
        f" sequence additions/removals"
    )
    if other:
        document.add_paragraph(
            f"{len(other)} further coding variant{'s' if len(other) != 1 else ''} with no "
            f"predicted protein change (e.g. synonymous), listed in the accompanying TSV"
        )
    document.add_paragraph(
        "All variants reported here passed quality filters (PASS), indicating high confidence."
    )

    # --- 2. Summary of findings ------------------------------------------
    document.add_heading("2. Summary of Findings", level=1)
    counts = []
    for sample in samples:
        per = Counter(v["category"] for v in gatk if v["sample"] == sample)
        sample_protein = [v for v in protein if v["sample"] == sample]
        label = str(per["protein_altering"])
        if sample_protein:
            label += f" ({gene_note(sample_protein).strip(' ()')})"
        counts.append([sample, label, per["regulatory"], per["indel"]])
    add_table(document,
              ["Sample", "Protein-Altering SNPs", "Regulatory SNPs", "Indels"],
              counts)

    # --- 3. Detailed SNP analysis ----------------------------------------
    document.add_heading("3. Detailed SNP Analysis", level=1)

    document.add_paragraph("A. Protein-Altering SNP (High Impact)", style="Subtitle")
    if protein:
        rows = [[v["sample"], v["pos"], f"{v['ref']} -> {v['alt']}", v["gene"],
                 v["annotation"].replace("_variant", "").replace("_", " ").capitalize(),
                 v["protein_change"]] for v in sorted(protein, key=lambda v: (v["pos"], v["sample"]))]
        add_table(document,
                  ["Sample", "Position", "Base change", "Gene", "Effect", "Amino Acid Change"],
                  rows)
    else:
        document.add_paragraph("No protein-altering SNPs passed filtering in these samples.")

    document.add_paragraph("B. Regulatory SNPs", style="Subtitle")
    if regulatory:
        rows = []
        for samples_joined, key in group_samples(
                sorted(regulatory, key=lambda v: (v["pos"], v["sample"])),
                ["pos", "ref", "alt", "gene", "annotation"]):
            pos, ref, alt, gene, annotation = key
            rows.append([samples_joined, pos, f"{ref} -> {alt}", gene,
                         annotation.replace("_gene_variant", "").capitalize()])
        add_table(document, ["Sample", "Position", "Base change", "Gene", "Region"], rows)
    else:
        document.add_paragraph("No regulatory SNPs passed filtering in these samples.")

    document.add_paragraph("C. Insertion/Deletion (Indel)", style="Subtitle")
    if indels:
        rows = []
        for v in sorted(indels, key=lambda v: (v["sample"], v["pos"])):
            rows.append([v["sample"], v["pos"],
                         f"{shorten(v['ref'])} -> {shorten(v['alt'])}",
                         indel_type(v["ref"], v["alt"]), v["gene"],
                         v["annotation"].replace("_gene_variant", "").capitalize()])
        add_table(document, ["Sample", "Position", "Changes", "Type", "Gene", "Region"], rows)
        document.add_paragraph(
            "Note: refer to the complete annotation details in the accompanying TSV file, "
            "as the full insertion sequences are too long to display in this summary table."
        )
    else:
        document.add_paragraph("No indels passed filtering in these samples.")

    # --- 4. Minor variants ------------------------------------------------
    document.add_heading("4. Key Summary - Minor variant analysis", level=1)
    document.add_paragraph(
        "The previous analysis (SNP) used strict quality filters to find high-confidence, "
        "common variants. In contrast, this analysis uses sensitive statistical models "
        "specifically designed to detect minor variants, rare mutations present at low "
        "frequencies in the viral population."
    )
    document.add_paragraph(
        f"Minor Variants Identified (allele frequency >= {args.min_af:.1%}, "
        f"{len(minor)} positions)"
    )
    rows = []
    for v in sorted(minor, key=lambda v: (v["sample"], v["pos"])):
        af = float(v["af"])
        depth = int(v["depth"]) if str(v["depth"]).isdigit() else v["depth"]
        rows.append([v["sample"], v["pos"], v["ref"], v["alt"],
                     f"{depth:,}" if isinstance(depth, int) else depth,
                     f"{af:.1%}"])
    add_table(document,
              ["Sample", "Position", "Reference", "Alternative", "Depth", "Allele Frequency"],
              rows,
              bold_predicate=lambda row: float(row[5].rstrip("%")) / 100 >= args.fixed_af)
    document.add_paragraph(
        f"Note: variants at >= {args.fixed_af:.1%} (bolded) are fixed mutations, not considered "
        f"as minor variants; these match the SNPs found in the earlier section. Complete data "
        f"is in the accompanying TSV file."
    )

    for paragraph in document.paragraphs:
        for run in paragraph.runs:
            if not run.font.size:
                run.font.size = Pt(11)
    return document


# ----------------------------------------------------------------------- main

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--proj-dir", required=True,
                        help="project directory, e.g. projects/SARS2-BavPat1")
    parser.add_argument("--samples", nargs="*", default=None,
                        help="sample names; default is every sample found")
    parser.add_argument("--out", default="variant_report",
                        help="output basename (.docx and .tsv are appended)")
    parser.add_argument("--reference", default="MZ558051.1 (SARS-CoV-2 BavPat1)",
                        help="reference name, printed in the report text")
    parser.add_argument("--title", default="SARS-CoV-2 variant report",
                        help="document title; pass '' for none")
    parser.add_argument("--min-af", type=float, default=0.02,
                        help="minimum allele frequency for the minor-variant table")
    parser.add_argument("--fixed-af", type=float, default=0.90,
                        help="allele frequency at or above which a variant is called fixed")
    parser.add_argument("--keep-nonpass", action="store_true",
                        help="include records that did not PASS filtering")
    args = parser.parse_args()

    variants_dir = os.path.join(args.proj_dir, "results", "02.variants")
    if not os.path.isdir(variants_dir):
        sys.exit(f"not a directory: {variants_dir}")

    samples = args.samples
    if not samples:
        found = glob.glob(os.path.join(variants_dir, "*.combined.annotate.vcf"))
        samples = sorted(os.path.basename(p).split(".combined.annotate.vcf")[0] for p in found)
    if not samples:
        sys.exit(f"no *.combined.annotate.vcf files in {variants_dir}")

    pass_only = not args.keep_nonpass
    all_variants = []
    for sample in samples:
        annotated = os.path.join(variants_dir, f"{sample}.combined.annotate.vcf")
        if not os.path.exists(annotated):
            sys.exit(f"missing {annotated}")
        gatk = read_annotated_vcf(annotated, sample, pass_only)
        minor_path = find_minor_vcf(variants_dir, sample)
        minor = read_minor_vcf(minor_path, sample, args.min_af, pass_only) if minor_path else []
        if minor_path is None:
            print(f"  ! no minor-variant VCF found for {sample}", file=sys.stderr)
        counts = Counter(v["category"] for v in gatk)
        print(f"{sample}: {len(gatk)} high-confidence "
              f"({counts['protein_altering']} protein-altering, {counts['regulatory']} regulatory, "
              f"{counts['indel']} indel, {counts['other']} other), {len(minor)} minor"
              + (f"   [{os.path.basename(minor_path)}]" if minor_path else ""))
        all_variants.extend(gatk)
        all_variants.extend(minor)

    columns = ["sample", "source", "chrom", "pos", "ref", "alt", "qual", "filter",
               "gene", "annotation", "impact", "protein_change", "category", "depth", "af"]
    tsv_path = f"{args.out}.tsv"
    with open(tsv_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for v in sorted(all_variants, key=lambda v: (v["sample"], v["source"], v["pos"])):
            writer.writerow(v)

    document = build_document(samples, all_variants, args)
    docx_path = f"{args.out}.docx"
    document.save(docx_path)
    print(f"\nwrote {docx_path} and {tsv_path}")


if __name__ == "__main__":
    main()
