## 003 - Trimming and QC of trimmed reads

from snakemake.io import glob_wildcards

# Discover paired-end samples
# SAMPLES = glob_wildcards(config["proj_dir"] + "/rawdata/{sample}_R1_001.fastq.gz").sample

# SAMPLES = glob_wildcards(config["data_dir"] + "/rawdata/P1658_{sample}_L002_R1_001_1.fq.gz").sample
# eg. P1658_001_JLEG2_12_S2_L002_R1_001_1.fq.gz

rule all003:
    input:
        expand(config["proj_dir"] + "/trimmed/{sample}_R1_tr_1P.fq.gz", sample=SAMPLES),
        expand(config["proj_dir"] + "/trimmed/{sample}_R2_tr_2P.fq.gz", sample=SAMPLES),
        config["proj_dir"] + "/QC_trimmed/final_multiqc_report.html"

rule trimmomatic:
    input: 
        R1 = RAW_R1,
        R2 = RAW_R2
    output:
        trim_1P = config["proj_dir"] + "/trimmed/{sample}_R1_tr_1P.fq.gz", 
        trim_1U = config["proj_dir"] + "/trimmed/{sample}_R1_tr_1U.fq.gz",  
        trim_2P = config["proj_dir"] + "/trimmed/{sample}_R2_tr_2P.fq.gz",  
        trim_2U = config["proj_dir"] + "/trimmed/{sample}_R2_tr_2U.fq.gz"   
    threads: config["threads"]
    params: config["adapter"]
    shell: 
        """
        mkdir -p {config[proj_dir]}/trimmed
        trimmomatic PE -threads {threads} \
        -Xmx16G \
        {input.R1} {input.R2} \
        {output.trim_1P} {output.trim_1U} \
        {output.trim_2P} {output.trim_2U} \
        ILLUMINACLIP:{params}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """

# FastQC on trimmed files
rule fastqc_trimmed:
    input: 
        config["proj_dir"] + "/trimmed/{sample}_{read}.fq.gz"
    output: 
        html = config["proj_dir"] + "/QC_trimmed/{sample}_{read}_fastqc.html",
        zip = config["proj_dir"] + "/QC_trimmed/{sample}_{read}_fastqc.zip"
    shell:
        """
        mkdir -p {config[proj_dir]}/QC_trimmed/
        fastqc {input} --outdir {config[proj_dir]}/QC_trimmed/
        """

# MultiQC on trimmed QC results
rule multiqc_trimmed:
    input:
        expand(config["proj_dir"] + "/QC_trimmed/{sample}_{read}_fastqc.html", 
               sample=SAMPLES, 
               read=["R1_tr_1P", "R2_tr_2P"])  # Only the paired outputs
    output:
        config["proj_dir"] + "/QC_trimmed/final_multiqc_report.html"
    params:
        report_name = "final_multiqc_report.html"     
    shell:
        "multiqc {config[proj_dir]}/QC_trimmed/ -o {config[proj_dir]}/QC_trimmed/ -n {params.report_name}"
