## To filter and annotate minor variants using LoFreq

from snakemake.io import glob_wildcards

# Discover paired-end samples
# SAMPLES = glob_wildcards(config["proj_dir"] + "/rawdata/{sample}_R1_001.fastq.gz").sample

# SAMPLES = glob_wildcards(config["data_dir"] + "/rawdata/P1658_{sample}_L002_R1_001_1.fq.gz").sample

rule all007:
    input:
        expand(config["proj_dir"] + "/results/02.variants/{sample}.lofreq.filtered.ann.vcf", sample=SAMPLES)

rule lofreq_indelqual:
    input:
        bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bam",
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        indelq_bam = config["proj_dir"] + "/results/02.variants/{sample}.lowfreq.indelq.bam",
        indelq_bai = config["proj_dir"] + "/results/02.variants/{sample}.lowfreq.indelq.bam.bai"
    shell:
        """
        lofreq indelqual --dindel -f {input.ref} -o {output.indelq_bam} {input.bam}
        samtools index {output.indelq_bam}
        """


rule lofreq_call:
    input:
        bam = rules.lofreq_indelqual.output.indelq_bam,
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.lofreq.vcf.gz",
        tbi = config["proj_dir"] + "/results/02.variants/{sample}.lofreq.vcf.gz.tbi"
    shell:
        """
        lofreq call -f {input.ref} -o {output.vcf} {input.bam}
        # overwrite existing .tbi if it exists
        tabix -f -p vcf {output.vcf}
        """


rule lofreq_filter:
    input:
        vcf = rules.lofreq_call.output.vcf
    output:
        filtered_vcf = config["proj_dir"] + "/results/02.variants/{sample}.minor_variants.vcf.gz",
        filtered_tbi = config["proj_dir"] + "/results/02.variants/{sample}.minor_variants.vcf.gz.tbi"
    params:
        tmp = config["proj_dir"] + "/results/02.variants/{sample}.lofreq_prefilter.tmp.vcf.gz"
    shell:
        """
        # first apply LoFreq filtering thresholds
        lofreq filter --af-min 0.005 --sb-alpha 0.05 --cov-min 100 \
            -i {input.vcf} -o {params.tmp}

        # then apply a stricter bcftools filter for coverage and AF cutoff
        bcftools view -i 'INFO/DP>=100 && INFO/AF>=0.03' \
            {params.tmp} -Oz -o {output.filtered_vcf}

        tabix -f -p vcf {output.filtered_vcf}
        rm -f {params.tmp} {params.tmp}.tbi
        """


rule bcftools_reheader:
    input:
        vcf = rules.lofreq_filter.output.filtered_vcf,
        fai = config["ref_dir"] + config["ref_fasta"] + ".fai"
    output:
        reheader_vcf = config["proj_dir"] + "/results/02.variants/{sample}.minor_variants.reheader.vcf.gz",
        reheader_tbi = config["proj_dir"] + "/results/02.variants/{sample}.minor_variants.reheader.vcf.gz.tbi"
    shell:
        """
        bcftools reheader --fai {input.fai} -o {output.reheader_vcf} {input.vcf}
        tabix -f -p vcf {output.reheader_vcf}
        """


rule annotate_lofreq:
    input:
        vcf = rules.bcftools_reheader.output.reheader_vcf
    output:
        ann_vcf = config["proj_dir"] + "/results/02.variants/{sample}.lofreq.filtered.ann.vcf",
        stats = config["proj_dir"] + "/results/02.variants/{sample}.snpEff_summary.html"
    params:
        snpjar = os.path.expanduser("~/apps/snpEff_custom/snpEff.jar"),
        genome = config["proj"]
    shell:
        """
        java -Xmx4g -jar {params.snpjar} -v -stats {output.stats} \
            {params.genome} {input.vcf} > {output.ann_vcf}
        """


### If neccessary, summarise the minor variants into a tsv with this commmand.
# bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/DP\t%INFO/AF\n' CHIKV-4.lofreq.filtered.ann.vcf > CHIKV-4.lofreq.filtered.summary.tsv
