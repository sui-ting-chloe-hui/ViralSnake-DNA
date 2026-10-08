# ViralSnake-DNA: from raw reads to annotated viral variants with Snakemake

A Snakemake pipeline that takes paired-end Illumina reads from a virus sample and finds the mutations in it, both the main (consensus) variants and low-frequency (minor) variants present in a small fraction of the viral population.

It has been used for CHIKV, SARS-CoV-2 and HIV-1 NL4-3. Any virus with a reference genome can be added: you put your files in place, write one config file, and run one command. No programming is needed.

---

## How it works

```
raw reads → QC → trimming → host filtering → align to virus → variant calling → annotation
                                                            └→ minor variants (LoFreq) → annotation
```

| Step | File | What it does |
|------|------|--------------|
| 001 | `001_qc.smk` | Checks the quality of the raw reads (FastQC) and summarises all samples in one report (MultiQC) |
| 002 | `002_index.smk` | Prepares the viral reference so reads can be aligned to it |
| 003 | `003_trim.smk` | Removes adapters and low-quality bases (Trimmomatic), then re-checks quality |
| 004a | `004a_host_filter.smk` | Removes reads from the host cells, depending on `host_mode` (see [Host filtering](#4-choose-how-to-handle-host-reads)) |
| 004b | `004b_align_virus.smk` | Aligns reads to the viral genome, marks PCR duplicates, and reports mapping statistics and coverage |
| 005b | `005b_variants.smk` | Calls variants |
| 006 | `006_annotateSNP.smk` | Annotates variants with their effect on genes and proteins (snpEff) |
| 007 | `007_minorvariants.smk` | Calls low-frequency variants with LoFreq, filters and annotates them |

Snakemake works out which steps need to run. If you stop halfway and restart, it carries on from where it stopped, and it never redoes a step whose output already exists and is up to date.

---

## Key terms

- **FASTQ**: the raw sequencing file. Paired-end sequencing gives two files per sample, **R1** and **R2**, one for each end of the DNA fragment.
- **Reference genome**: a FASTA file with the known sequence of your virus. Reads are compared against it to find differences.
- **Alignment (BAM file)**: where each read belongs on the reference.
- **Host reads**: reads from the cells the virus was grown in (or from the patient), not from the virus.
- **Variant (VCF file)**: a position where your sample differs from the reference.
- **Consensus variant**: a change carried by most of the viral population.
- **Minor variant**: a change carried by only a small fraction, e.g. 3–20% of reads.
- **AF (allele frequency)**: the fraction of reads carrying the variant.
- **DP (depth)**: how many reads cover that position. More depth means more confidence.
- **Annotation**: what a variant does, e.g. synonymous, missense (amino acid change) or stop codon.

---

## 1. Install the software

The simplest way is one conda environment:

```
conda create -n viral-variants -c conda-forge -c bioconda \
    snakemake fastqc multiqc trimmomatic bwa samtools bcftools htslib lofreq gatk4 openjdk
conda activate viral-variants
```

If conda struggles to solve the environment, create one environment for Snakemake and a second one for the tools, and activate both (`conda activate --stack`).

**snpEff** is used with a custom database for each virus. The pipeline expects it at `~/apps/snpEff_custom/snpEff.jar`. If yours is elsewhere, change the `snpjar` path in `workflows/006_annotateSNP.smk` and `workflows/007_minorvariants.smk`.

---

## 2. Get the pipeline

```
git clone https://github.com/YOUR_USERNAME/viral-variant-pipeline.git
cd viral-variant-pipeline
```

**Always run Snakemake from this top folder.** All paths in the config files are relative to it.

---

## 3. Set up your project

Each virus or experiment is a "project". This example uses the name `MYVIRUS`; replace it with your own (no spaces).

### 3.1 Create the folders

```
mkdir -p projects/MYVIRUS/rawdata projects/MYVIRUS/reference
```

`projects/` is never uploaded to GitHub, so your data stays on your server.

### 3.2 Add your raw reads

Copy or link your FASTQ files into `projects/MYVIRUS/rawdata/`. Linking saves space:

```
ln -s /path/to/sequencing/run/*.fastq.gz /path/to/sequencing/run/*.fastq.gz.md5 projects/MYVIRUS/rawdata/

# Optional step to check for the integrity of the FASTQ files:
cd projects/SARS2-BavPat1/rawdata
for f in *.md5; do md5sum -c "$f"; done

# Every line should end in OK. If one says FAILED, download that file again.

```

You can also leave the files where they are and point `raw_dir` at them in the config (step 3.5).

### 3.3 Add the reference genome

Put the viral reference **FASTA** (sequence) and **GFF** (gene annotation) in `projects/MYVIRUS/reference/`. Both are usually available from NCBI for the strain or plasmid you used.

Make the sequence dictionary once (step 002 creates the other index files):

```seqki
samtools dict projects/MYVIRUS/reference/myvirus.fasta \
    -o projects/MYVIRUS/reference/myvirus.fasta.dict
```

The name must be `<fasta file name>.dict`, e.g. `myvirus.fasta.dict`, not `myvirus.dict`.

### 3.4 Build a snpEff database

snpEff needs a database built from your FASTA and GFF. For downloading the GFF, the easiest route is NCBI's sequence viewer CGI, which will render any annotated nuccore record as GFF3. Then use 'curl' to retrieve the file.

 **The database name must be the same as your project name** (`MYVIRUS`). You can follow the snpEff documentation on building a database for a custom genome. Here's what I do for building the database:

```
Steps (one-time setup):
#   1. Locate your snpEff installation:
        find / -type f -name "snpEff.jar" 2>/dev/null
#   2. Copy it to your user directory:
         mkdir -p ~/apps
         cp -r /opt/miniconda/envs/ntm-profiler/share/snpeff-5.1-1 ~/apps/snpEff_custom
#   3. Edit ~/apps/snpEff_custom/snpEff.config (use nano or vi):
#      Add a line for your genome under 'genomes':
#         MYVIRUS.genome : MYVIRUS_reference
#   4. Place genome files in the data directory of snpEff:
         ~/apps/snpEff_custom/data/MYVIRUS/sequences.fa
         ~/apps/snpEff_custom/data/MYVIRUS/genes.gff
#   5. Build the database once:
        cd ~/apps/snpEff_custom
        java -jar snpEff.jar build -gff3 -noGenome -noCheckCds -noCheckProtein -v MYVIRUS
# p.s. (-gtf22) if using GTF format 
```

### 3.5 Write the config file

Copy an existing config and edit it:

```
cp config/CHIKV.config.yaml config/MYVIRUS.config.yaml
```

```yaml
proj: "MYVIRUS"
proj_dir: "projects/MYVIRUS"
ref_dir: "projects/MYVIRUS/reference/"
ref_fasta: "myvirus.fasta"

raw_dir: "projects/MYVIRUS/rawdata"
r1_pattern: "{sample}_R1_001.fastq.gz"
r2_pattern: "{sample}_R2_001.fastq.gz"

adapter: "/path/to/trimmomatic/adapters/TruSeq3-PE.fa"
threads: 8

host_mode: "remove"
host_ref: "Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
host_ref_dir: "/path/to/reference/hg38/"
min_mapq: 20
```

| Setting | What to put |
|---------|-------------|
| `proj` | Your project name. Must match the config file name and the snpEff database name. |
| `proj_dir` | Project folder. All outputs go here. No slash at the end. |
| `ref_dir` | Folder with the reference files. **Keep the slash at the end.** |
| `ref_fasta` | Name of the reference FASTA file inside `ref_dir`. |
| `raw_dir` | Folder with the raw FASTQ files. Can be anywhere on the server. |
| `r1_pattern`, `r2_pattern` | How your read files are named, with `{sample}` where the sample name goes (see below). |
| `adapter` | Trimmomatic adapter file. Its path changes when Trimmomatic is updated, so check it if trimming fails. |
| `threads` | CPU threads per job. Snakemake never uses more than the `--cores` you give it. |
| `host_mode` | `none`, `remove` or `competitive`. See section 4. |
| `host_ref`, `host_ref_dir` | Host genome FASTA and its folder. Not needed if `host_mode` is `none`. |
| `min_mapq` | Mapping-quality cut-off, used only by `competitive` mode. |

**Matching your file names.** `{sample}` marks the part of the file name that differs between samples. Everything else must be identical across files.

| File | Pattern | Sample name found |
|------|---------|-------------------|
| `CHIKV-1_R1_001.fastq.gz` | `{sample}_R1_001.fastq.gz` | `CHIKV-1` |
| `P1658_001_JLEG2_1_S1_L002_R1_001_1.fq.gz` | `P1658_{sample}_L002_R1_001_1.fq.gz` | `001_JLEG2_1_S1` |

The sample name is used to name every output file.

---

## 4. Choose how to handle host reads

Your samples contain DNA from the cells the virus grew in (or from the patient). Most host reads simply fail to align to the small viral genome, but a few look enough like the virus to be forced onto it. They then appear as false low-frequency variants, which is exactly what the minor-variant step looks for.

| `host_mode` | What happens | When to use it |
|-------------|--------------|----------------|
| `none` | No filtering. Trimmed reads go straight to viral alignment. | Samples with very little host DNA (e.g. plasmid or purified virus), or when you only need consensus variants. |
| `remove` | Reads are aligned to the host genome first; pairs that match the host are removed. | **Recommended default.** Reliable and fast after the host index is built once. |
| `competitive` | Reads are aligned to host + virus together; only reads that clearly belong to the virus are kept. | When a virus resembles its host's genome (e.g. HIV and human endogenous retroviruses). Most thorough, but builds a large index for each virus. |

In every mode, the kept reads are then aligned to the viral genome only, so the variant steps work the same way.

**The host does not have to be human.** Set `host_ref` to the genome of the cells you used, for example *Chlorocebus sabaeus* for Vero cells, or an *Aedes* genome for mosquito cell lines.

**The first run is slow.** The host genome has to be indexed (for a human-sized genome, roughly 1–2 hours and ~5 GB of memory). The index is reused afterwards. In `competitive` mode, a combined index is built for each virus.

**HIV and `competitive` mode.** HIV's two LTRs are identical, so reads there cannot be placed uniquely and get a low mapping quality. With `min_mapq: 20`, they are removed. Use `remove` mode, or set `min_mapq: 0`, if you need LTR coverage.

---

## 5. Run the pipeline

### 5.1 Dry run first (always)

A dry run shows what would happen without running anything:

```
snakemake -s workflows/Snakefile_main --config proj=MYVIRUS -n
```

At the bottom you'll see a table of jobs. For a new project, check that the counts make sense: e.g. with 10 samples you should see about 20 `raw_fastqc` jobs (R1 and R2 for each). If it says **"Nothing to be done"** for a brand-new project, no samples were found, so check `raw_dir` and the file patterns.

### 5.2 Real run

Runs can take hours, so start them inside `tmux` (or `screen`). That way they keep going if you close your laptop or lose connection:

```
tmux new -s myvirus
conda activate viral-variants
snakemake -s workflows/Snakefile_main --config proj=MYVIRUS --cores 8 -p
```

Detach with `Ctrl+B` then `D`. Reattach later with `tmux attach -t myvirus`.

### 5.3 Run only part of the pipeline

Each step has its own target rule. For example, to run only QC:

```
snakemake -s workflows/Snakefile_main --config proj=MYVIRUS --cores 8 all001
```

Snakemake also runs any earlier steps a target depends on. To list the targets:

```
grep -n "^rule all" workflows/*.smk
```

### 5.4 Useful options

| Option | What it does |
|--------|--------------|
| `-n` | Dry run: show the plan, run nothing |
| `-p` | Print each command as it runs |
| `-r` | Show why each job is going to run |
| `--cores 8` | Maximum number of CPU cores to use |
| `--forceall -n` | List every job in the pipeline, including finished ones |
| `--rerun-triggers mtime` | Rerun only when input files are newer, not when code has been edited |
| `--unlock` | Remove the lock left by a run that was interrupted |

---

## 6. Check your results

Outputs are in `projects/MYVIRUS/`. Look at them in this order:

1. **Read quality**: `QC/raw_multiqc_report.html` and the report in `QC_trimmed/`. Open them in a browser. Look for low quality scores or leftover adapters.
2. **Mapping**: `results/01.alignment/<sample>.sorted.dedup.RG.flagstat.txt`. The "mapped" percentage tells you how much of each sample is virus.
3. **Coverage**: `results/01.alignment/<sample>.coveragedepth.txt`. Regions with low depth give unreliable variant calls.
4. **Variants**: the VCF files and per-sample snpEff reports (`<sample>.snpEff_summary.html`) under `results/`.
5. **Minor variants**: `results/02.variants/<sample>.lofreq.filtered.ann.vcf`.

| Folder | Contents |
|--------|----------|
| `QC/`, `QC_trimmed/` | Quality reports before and after trimming |
| `trimmed/` | Trimmed reads |
| `clean/` | Reads after host filtering |
| `results/01.alignment/` | Alignments (BAM), mapping statistics, coverage |
| `results/02.variants/` | Variant calls, minor variants and snpEff reports |
| `logs/` | Log files |

To turn a minor-variant VCF into a table you can open in Excel:

```
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n' \
    projects/MYVIRUS/results/02.variants/SAMPLE.lofreq.filtered.ann.vcf > SAMPLE.minor_variants.tsv
```

---

## 7. Build the variant report (optional)

`workflows/manual/variant_report.py` turns the variant calls into a Word
document with summary tables, plus a TSV holding every variant in one tidy
table. It replaces the `bcftools query` one-liner in section 6 when you want
something to hand to a collaborator.

It needs `python-docx`, which is not in the conda environment above:

```
pip install python-docx
```

(The PyPI package is `python-docx` but imports as `docx`. Installing `docx`
instead gets you an unrelated, abandoned package.)

```
python workflows/manual/variant_report.py
--proj-dir projects/MYVIRUS
--samples SAMPLE1 SAMPLE2
--reference "NC_000000.1 (MYVIRUS strain)"
--title "MYVIRUS variant report"
--out myvirus.report
```

This writes `myvirus.report.docx` and `myvirus.report.tsv`. Omit `--samples`
and every sample with a `combined.annotate.vcf` is picked up. `--reference` and
`--title` are printed into the report text only and default to SARS-CoV-2, so
set them for your own virus. Run `--help` for the remaining options, which
adjust the allele-frequency thresholds for the minor-variant table.


## 8. Troubleshooting

| Message or symptom | Likely cause and fix |
|--------------------|----------------------|
| `KeyError: 'something'` | A setting is missing from your config file. Compare it with an existing config. |
| `MissingInputException` | A file Snakemake needs doesn't exist. Read the path in the message: it is often the `.dict` file (step 3.3) or a typo in `ref_dir` / `ref_fasta`. |
| "Nothing to be done" on a new project | No samples found. Check `raw_dir`, `r1_pattern` and `r2_pattern` against `ls` of your raw data folder. |
| `LockException` | A previous run was interrupted. Run the same command with `--unlock`, then run again. |
| Trimming fails | The `adapter` path is probably out of date. Find the current one with `find ~/miniconda3 -name "TruSeq3-PE.fa"`. |
| Everything reruns unexpectedly | Code or settings changed. Check with `-n -r`; add `--rerun-triggers mtime` if the change shouldn't affect results. |
| `No such file or directory` for config | You are not in the pipeline's top folder, or the config name doesn't match `proj`. |

---

## Notes and limitations

- **Amplicon data.** Step 004b marks PCR duplicates. This is correct for shotgun libraries, but for amplicon (tiled PCR) sequencing, all reads in an amplicon start at the same position and would be treated as duplicates, removing most of the coverage. Do not use this pipeline unchanged for amplicon data.
- **Minor-variant thresholds** are set in the `lofreq_filter` rule of `007_minorvariants.smk`: LoFreq filtering (AF ≥ 0.005, strand-bias alpha 0.05, depth ≥ 100), then bcftools (depth ≥ 100, AF ≥ 0.03). Adjust them for your depth and question.
- **Contig names in `competitive` mode** must not be identical between the virus and host FASTA files.

---

## Repository layout

```
.
├── README.md
├── .gitignore
├── workflows/
│   ├── Snakefile_main        entry point; loads the config and all steps
│   ├── 001_qc.smk ... 007_minorvariants.smk
│   └── manual/               optional extras run by hand (e.g. coverage.smk)
├── config/
│   └── <PROJECT>.config.yaml one config per project
├── docs/                     workflow diagrams and notes
└── projects/                 your data and results (not uploaded to GitHub)
```
