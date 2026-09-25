## 004 - Alignment

from snakemake.io import glob_wildcards

# Discover paired-end samples from rawdata
# SAMPLES = glob_wildcards(config["proj_dir"] + "/rawdata/{sample}_R1_001.fastq.gz").sample
# SAMPLES = glob_wildcards(config["data_dir"] + "/rawdata/P1658_{sample}_L002_R1_001_1.fq.gz").sample

def reads_for_alignment(wildcards):
    if HOST_MODE == "none":
        r1, r2 = TRIM_R1, TRIM_R2
    else:
        r1, r2 = CLEAN_R1, CLEAN_R2
    return {"clean_R1": r1.format(sample=wildcards.sample),
            "clean_R2": r2.format(sample=wildcards.sample)}

rule all004b:
    input:
        expand(config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bam", sample=SAMPLES),
        expand(config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bai", sample=SAMPLES),
        expand(config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.flagstat.txt", sample=SAMPLES),
        expand(config["proj_dir"] + "/results/01.alignment/{sample}.coveragedepth.txt", sample=SAMPLES)
        # expand(config["proj_dir"] + "/results/01.alignment/{sample}.unmapped.fastq", sample=SAMPLES)

# Map reads to viral reference
rule align_to_reference:
    input:
        unpack(reads_for_alignment),
        ref_path = config["ref_dir"] + config["ref_fasta"]
    output:
        sorted_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.bam"
    threads: config["threads"]
    shell:
        """
        mkdir -p {config[proj_dir]}/results/01.alignment
        bwa mem -t {threads} -M {input.ref_path} {input.clean_R1} {input.clean_R2} | \
        samtools view -bS - | samtools sort -o {output.sorted_bam}
        samtools index {output.sorted_bam}
        """

# Mark duplicates
rule picard_mark_duplicates:
    input:
        sorted_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.bam"
    output:
        dedup_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.bam",
        metrics = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.metrics.txt"
    shell:
        """
        PicardCommandLine MarkDuplicates \
            I={input.sorted_bam} \
            O={output.dedup_bam} \
            M={output.metrics} \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT
        """

rule add_readgroups:
    input:
        dedup_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.bam"
    output:
        rg_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bam",
        rg_bai = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bai"
    shell:
        """
        gatk AddOrReplaceReadGroups \
            --INPUT {input.dedup_bam} \
            --OUTPUT {output.rg_bam} \
            --CREATE_INDEX true \
            --RGID {wildcards.sample} \
            --RGLB lib1 \
            --RGPL ILLUMINA \
            --RGPU unit1 \
            --RGSM {wildcards.sample}
        """

# Get alignment statistics (flagstat)
rule flagstat:
    input:
        final_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bam"
    output:
        stats = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.flagstat.txt"
    shell:
        """
        samtools flagstat {input.final_bam} > {output.stats}
        """

# Get depth of coverage information
rule depth_of_coverage:
    input:
        final_bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bam",
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        depth = config["proj_dir"] + "/results/01.alignment/{sample}.coveragedepth.txt"
    shell:
        """
        PicardCommandLine CollectWgsMetrics \
        -I {input.final_bam}\
        -O {output.depth}\
        -R {input.ref}
        """

# Extract unmapped reads (optional)
# rule extract_unmapped:
   # input:
    #    bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.bam"
   # output:
    #    unmapped_bam = config["proj_dir"] + "/results/01.alignment/{sample}.unmapped.bam",
     #   unmapped_fastq = config["proj_dir"] + "/results/01.alignment/{sample}.unmapped.fastq"
   # shell:
    #    """
     #   samtools view -b -f 4 {input.bam} > {output.unmapped_bam}
      #  samtools index {output.unmapped_bam}
       # samtools fastq {output.unmapped_bam} > {output.unmapped_fastq}
       # """



