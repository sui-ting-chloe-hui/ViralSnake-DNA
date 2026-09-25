#########################################
## 005 - Variant Calling with GATK
#########################################

from snakemake.io import glob_wildcards

# Discover paired-end samples
# SAMPLES = glob_wildcards(config["proj_dir"] + "/rawdata/{sample}_R1_001.fastq.gz").sample
#SAMPLES = glob_wildcards(config["data_dir"] + "/rawdata/P1658_{sample}_L002_R1_001_1.fq.gz").sample


rule all005b:
    input:
        expand(config["proj_dir"] + "/results/02.variants/{sample}.final_snps.vcf", sample=SAMPLES),
        expand(config["proj_dir"] + "/results/02.variants/{sample}.final_indels.vcf", sample=SAMPLES)


# GATK HaplotypeCaller
rule gatk_haplotypecaller:
    input: 
        bam = config["proj_dir"] + "/results/01.alignment/{sample}.sorted.dedup.RG.bam",
        ref = config["ref_dir"] + config["ref_fasta"]
    output: 
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw.vcf.gz"
    threads: config["threads"]
    shell:
        """
        mkdir -p {config[proj_dir]}/results/02.variants
        gatk --java-options "-Xmx16G" HaplotypeCaller \
            -R {input.ref} \
            -I {input.bam} \
            -O {output.vcf} \
            --pcr-indel-model NONE \
            --standard-min-confidence-threshold-for-calling 30 \
            --min-base-quality-score 20 \
            --sample-ploidy 1 \
            --annotation QualByDepth 
        """

rule extract_snps_quality:
    input:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw.vcf.gz",
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        snp_vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw_snps_recal.vcf.gz"
    threads: 1
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include SNP -O {output.snp_vcf}
        """

# Extract INDELs
rule extract_indels:
    input:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw.vcf.gz",
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        indel_vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw_indels_recal.vcf.gz"
    threads: 1
    shell:
        """
        gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-include INDEL -O {output.indel_vcf}
        """

# GATK SNP Filtering 
rule gatk_filter_snps:
    input:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw_snps_recal.vcf.gz",
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        filtered_snp = config["proj_dir"] + "/results/02.variants/{sample}.filtered_snps.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.filtered_snp} \
            -filter-name "QD_filter" -filter "QD < 2.0" \
	        -filter-name "FS_filter" -filter "FS > 60.0" \
	        -filter-name "MQ_filter" -filter "MQ < 40.0" \
	        -filter-name "SOR_filter" -filter "SOR > 4.0" \
	        -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	        -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	        -genotype-filter-expression "DP < 10" \
	        -genotype-filter-name "DP_filter" \
	        -genotype-filter-expression "GQ < 10" \
	        -genotype-filter-name "GQ_filter"
        """

rule gatk_variants_indel_filter:
    input:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.raw_indels_recal.vcf.gz",
        ref = config["ref_dir"] + config["ref_fasta"]
    output:
        filtered_indel = config["proj_dir"] + "/results/02.variants/{sample}.filtered_indel.vcf.gz"
    shell:
        """
        gatk VariantFiltration \
            -R {input.ref} \
            -V {input.vcf} \
            -O {output.filtered_indel} \
            -filter-name "QD_filter" -filter "QD < 2.0" \
	        -filter-name "FS_filter" -filter "FS > 200.0" \
	        -filter-name "SOR_filter" -filter "SOR > 10.0" \
	        -genotype-filter-expression "DP < 10" \
	        -genotype-filter-name "DP_filter" \
	        -genotype-filter-expression "GQ < 10" \
	        -genotype-filter-name "GQ_filter"
        """

# Create final SNPs VCF - that keep only the select variants that has PASS filter
rule create_final_snps:
    input:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.filtered_snps.vcf.gz"
    output:
        tmp_snpvcf = config["proj_dir"] + "/results/02.variants/{sample}.tmp_snps.vcf",
        final_snpvcf = config["proj_dir"] + "/results/02.variants/{sample}.final_snps.vcf"
    threads: 1
    shell:
        """
        # Select only PASS variants
        gatk SelectVariants --exclude-filtered -V {input.vcf} -O {output.tmp_snpvcf}
        
        # Remove genotype-level filter patterns
        cat {output.tmp_snpvcf} | grep -v -E "DP_filter|GQ_filter" > {output.final_snpvcf}
        """

# Create final INDELs VCF - allow multiple nucleotides
rule create_final_indels:
    input:
        vcf = config["proj_dir"] + "/results/02.variants/{sample}.filtered_indel.vcf.gz"
    output:
        tmp_indelvcf = config["proj_dir"] + "/results/02.variants/{sample}.tmp_indels.vcf",
        final_indelvcf = config["proj_dir"] + "/results/02.variants/{sample}.final_indels.vcf"
    threads: 1
    shell:
        """
        # Select only PASS variants
        gatk SelectVariants --exclude-filtered -V {input.vcf} -O {output.tmp_indelvcf}

        # Remove genotype-level filter patterns
        cat {output.tmp_indelvcf} | grep -v -E "DP_filter|GQ_filter" > {output.final_indelvcf}
        """