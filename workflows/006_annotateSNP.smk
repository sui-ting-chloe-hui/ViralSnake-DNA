#########################################
# 006 - Annotate SNPs using SnpEff
#########################################
# Before using SnpEff, make sure you have built your custom genome database.
# Steps (one-time setup):
#   1. Locate your snpEff installation:
#         find / -type f -name "snpEff.jar" 2>/dev/null
#   2. Copy it to your user directory:
#         mkdir -p ~/apps
#         cp -r /opt/miniconda/envs/ntm-profiler/share/snpeff-5.1-1 ~/apps/snpEff_custom
#   3. Edit ~/apps/snpEff_custom/snpEff.config (use nano or vi):
#      Add a line for your genome under 'genomes' (example here is SARS2, replace with your own genome):
#         SARS2.genome : SARS-CoV-2_reference
#   4. Place genome files in the data directory of snpEff:
#         ~/apps/snpEff_custom/data/SARS2/sequences.fa
#         ~/apps/snpEff_custom/data/SARS2/genes.gff
#   5. Build the database once:
#         cd ~/apps/snpEff_custom
#         java -jar snpEff.jar build -gff3 -noGenome -noCheckCds -noCheckProtein -v SARS2
# p.s. (-gtf22) if using GTF format 
#########################################

from snakemake.io import glob_wildcards

# Discover paired-end samples
# SAMPLES = glob_wildcards(config["proj_dir"] + "/rawdata/{sample}_R1_001.fastq.gz").sample
# SAMPLES = glob_wildcards(config["data_dir"] + "/rawdata/P1658_{sample}_L002_R1_001_1.fq.gz").sample

rule all006:
    input:
        expand(config["proj_dir"] + "/results/02.variants/{sample}.annotate_snps.vcf", sample=SAMPLES),
        expand(config["proj_dir"] + "/results/02.variants/{sample}.annotate_indels.vcf", sample=SAMPLES),
        expand(config["proj_dir"] + "/results/03.consensus/{sample}.consensus.fasta", sample=SAMPLES)

rule snpEff_annotate_snps:
    input:
        filtered_vcf = config["proj_dir"] + "/results/02.variants/{sample}.final_snps.vcf"
    output:
        annotate_vcf = config["proj_dir"] + "/results/02.variants/{sample}.annotate_snps.vcf"
    params:
        snpeff_jar = "~/apps/snpEff_custom/snpEff.jar",
        snpeff_config = "~/apps/snpEff_custom/snpEff.config",
        genome = config["proj"] ##SARS2-3 needs to change genome to "SARS2GFP"
    log:
        config["proj_dir"] + "/logs/snpEff/{sample}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        cat {input.filtered_vcf} | java -Xmx4g -jar {params.snpeff_jar} \
            -v {params.genome} \
            -c {params.snpeff_config} \
            - > {output.annotate_vcf} 2> {log}
        """

## sSNPs
#for file in SARS2-*.annotate_snps.vcf; do
#grep 'synonymous_variant' "$file" | awk -F'\t' '{split($5, alt, ","); split($8, ann, "|"); print $1, $2, $4, alt[1], ann[3], ann[4], ann[11], ann[12], ann[13], ann[14], ann[15], $7}' > "${file%.annotate_snps.vcf}.snpeff.sS.txt"
#done

## nsSNPs
#for file in SARS2-*.annotate_snps.vcf; do
#grep 'missense_variant' "$file" | awk -F'\t' '{split($5, alt, ","); split($8, ann, "|"); print $1, $2, $4, alt[1], ann[3], ann[4], ann[11], ann[12], ann[13], ann[14], ann[15], $7}' > "${file%.annotate_snps.vcf}.snpeff.nS.txt"
#done

#for file in SARS2-*.annotate_indels.vcf; do
#grep -E 'frameshift_variant|inframe_insertion|inframe_deletion' "$file" | awk -F'\t' '{split($5, alt, ","); split($8, ann, "|"); print $1, $2, $4, alt[1], ann[3], ann[4], ann[11], ann[12], ann[13], ann[14], ann[15], $7}' > "${file%.annotate_indels.vcf}.snpeff.indels.txt"
#done


rule snpEff_annotate_indels:
    input:
        filtered_vcf = config["proj_dir"] + "/results/02.variants/{sample}.final_indels.vcf"
    output:
        annotate_vcf = config["proj_dir"] + "/results/02.variants/{sample}.annotate_indels.vcf"
    params:
        snpeff_jar = "~/apps/snpEff_custom/snpEff.jar",
        snpeff_config = "~/apps/snpEff_custom/snpEff.config",
        genome = config["proj"] ##SARS2-3 needs to change genome to SARS2GFP
    log:
        config["proj_dir"] + "/logs/snpEff/{sample}.log"
    shell:
        """
        mkdir -p $(dirname {log})
        cat {input.filtered_vcf} | java -Xmx4g -jar {params.snpeff_jar} \
            -v {params.genome} \
            -c {params.snpeff_config} \
            - > {output.annotate_vcf} 2> {log}
        """

rule cat_vcf:
    input:
        snp_vcf = config["proj_dir"] + "/results/02.variants/{sample}.annotate_snps.vcf",
        indel_vcf = config["proj_dir"] + "/results/02.variants/{sample}.annotate_indels.vcf"
    output:
        combined_vcf = config["proj_dir"] + "/results/02.variants/{sample}.combined.annotate.vcf"
    shell:
        """
        gatk MergeVcfs -I {input.snp_vcf} -I {input.indel_vcf} -O {output.combined_vcf}
        """

rule index_combined_vcf:
    input:
        combined_vcf = config["proj_dir"] + "/results/02.variants/{sample}.combined.annotate.vcf"
    output:
        index = config["proj_dir"] + "/results/02.variants/{sample}.combined.annotate.vcf.idx"
    shell:
        """
        gatk IndexFeatureFile -I {input.combined_vcf}
        """

rule generate_consensus_sequence:
    input:
        combined_vcf = config["proj_dir"] + "/results/02.variants/{sample}.combined.annotate.vcf",
        ref_fasta = config["ref_dir"] + config["ref_fasta"],
        ref_dict = (config["ref_dir"] + config["ref_fasta"]).rsplit(".", 1)[0] + ".dict",
        ref_fai = config["ref_dir"] + config["ref_fasta"] + ".fai"
    output:
        consensus_fasta = config["proj_dir"] + "/results/03.consensus/{sample}.consensus.fasta"
    shell:
        """
        mkdir -p {config[proj_dir]}/results/03.consensus
        gatk FastaAlternateReferenceMaker -R {input.ref_fasta} -V {input.combined_vcf} -O {output.consensus_fasta}
        """
## consenus sequence
# bcftools annotate -x FILTER -O v -o 

# gatk IndexFeatureFile CHIKV-1.annotate.vcf
# gatk FastaAlternateReferenceMaker -R ../../reference/CHIKV_NC2011_568.fasta -V CHIKV-1.annotate.vcf -O CHIKV1.consensus.fasta
