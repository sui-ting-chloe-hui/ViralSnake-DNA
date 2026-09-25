## 004a - Host read filtering
##
## host_mode in the config decides what happens before viral alignment (004b):
##   "none"        (A) no host filtering; 004b aligns the trimmed reads directly
##   "remove"      (B) align to the host genome; keep read pairs where both mates are unmapped
##   "competitive" (C) align to a combined host + virus reference; keep read pairs where
##                     both mates map to the virus with MAPQ >= min_mapq
##
## B and C write the filtered reads to clean/. 004b then aligns them to the viral
## reference only, so every downstream step (variants, snpEff, LoFreq) is unchanged.
##
## The host can be any organism: set host_ref / host_ref_dir in the config
## (e.g. human, Vero / African green monkey, mosquito cell lines).


# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

HOST_MODE = config.get("host_mode", "remove")
if HOST_MODE not in ("none", "remove", "competitive"):
    raise ValueError("host_mode must be 'none', 'remove' or 'competitive', got: " + str(HOST_MODE))

# host_ref / host_ref_dir; falls back to the older human_ref / human_ref_dir keys
HOST_FASTA = (config.get("host_ref_dir", config.get("human_ref_dir", ""))
              + config.get("host_ref", config.get("human_ref", "")))

VIRAL_FASTA = config["ref_dir"] + config["ref_fasta"]
COMPETITIVE_DIR = config["proj_dir"] + "/reference/competitive"
COMBINED_FASTA = COMPETITIVE_DIR + "/host_plus_" + config["proj"] + ".fa"
VIRAL_CONTIGS = COMPETITIVE_DIR + "/viral_contigs.txt"
MIN_MAPQ = config.get("min_mapq", 20)

BWA_EXT = ["amb", "ann", "bwt", "pac", "sa"]

# Read files going in (from 003) and coming out (used by 004b)
TRIM_R1 = config["proj_dir"] + "/trimmed/{sample}_R1_tr_1P.fq.gz"
TRIM_R2 = config["proj_dir"] + "/trimmed/{sample}_R2_tr_2P.fq.gz"
CLEAN_R1 = config["proj_dir"] + "/clean/{sample}_R1_clean.fastq.gz"
CLEAN_R2 = config["proj_dir"] + "/clean/{sample}_R2_clean.fastq.gz"


rule all004a:
    input:
        [] if HOST_MODE == "none"
        else expand(CLEAN_R1, sample=SAMPLES) + expand(CLEAN_R2, sample=SAMPLES)


# ---------------------------------------------------------------------------
# B: remove reads that map to the host
# ---------------------------------------------------------------------------

if HOST_MODE == "remove":

    rule bwa_index_host:
        input:
            HOST_FASTA
        output:
            expand(HOST_FASTA + ".{ext}", ext=BWA_EXT)
        shell:
            "bwa index {input}"

    rule host_remove:
        input:
            r1 = TRIM_R1,
            r2 = TRIM_R2,
            ref = HOST_FASTA,
            idx = expand(HOST_FASTA + ".{ext}", ext=BWA_EXT)
        output:
            r1 = CLEAN_R1,
            r2 = CLEAN_R2
        threads: config.get("threads", 8)
        shell:
            """
            # -f 12: both mates unmapped to the host; -F 0x900: no secondary/supplementary
            bwa mem -t {threads} {input.ref} {input.r1} {input.r2} \
                | samtools fastq -n -f 12 -F 0x900 \
                    -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -
            """


# ---------------------------------------------------------------------------
# C: competitive mapping to host + virus
# ---------------------------------------------------------------------------

if HOST_MODE == "competitive":

    rule combined_reference:
        input:
            host = HOST_FASTA,
            virus = VIRAL_FASTA
        output:
            COMBINED_FASTA
        shell:
            """
            if [[ {input.host} == *.gz ]]; then zcat {input.host}; else cat {input.host}; fi > {output}
            cat {input.virus} >> {output}
            """

    rule viral_contig_list:
        input:
            VIRAL_FASTA
        output:
            VIRAL_CONTIGS
        shell:
            "grep '^>' {input} | sed 's/^>//; s/[[:space:]].*//' > {output}"

    rule bwa_index_combined:
        input:
            COMBINED_FASTA
        output:
            expand(COMBINED_FASTA + ".{ext}", ext=BWA_EXT)
        shell:
            "bwa index {input}"

    rule competitive_filter:
        input:
            r1 = TRIM_R1,
            r2 = TRIM_R2,
            ref = COMBINED_FASTA,
            idx = expand(COMBINED_FASTA + ".{ext}", ext=BWA_EXT),
            contigs = VIRAL_CONTIGS
        output:
            r1 = CLEAN_R1,
            r2 = CLEAN_R2
        params:
            mapq = MIN_MAPQ
        threads: config.get("threads", 8)
        shell:
            """
            # keep primary alignments on viral contigs with MAPQ >= min_mapq,
            # then keep only pairs where both mates passed
            bwa mem -t {threads} {input.ref} {input.r1} {input.r2} \
                | samtools view -h -q {params.mapq} -F 0x904 - \
                | awk 'BEGIN {{ while ((getline c < "{input.contigs}") > 0) viral[c] = 1 }} /^@/ || ($3 in viral)' \
                | samtools collate -O -u - \
                | samtools fastq -n -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -
            """
