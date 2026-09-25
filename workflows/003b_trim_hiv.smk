from snakemake.io import glob_wildcards

# Capture everything between P1658_ and _L002_R1...
# SAMPLES_FULL, = glob_wildcards(config["data_dir"] + "/rawdata/P1658_{sample_full}_L002_R1_001_1.fq.gz")

rule all003b:
    input:
        expand(config["proj_dir"] + "/trimmed/{sample}_R1_tr_1P.fq.gz", sample=SAMPLES_FULL),
        expand(config["proj_dir"] + "/trimmed/{sample}_R2_tr_2P.fq.gz", sample=SAMPLES_FULL),
        config["proj_dir"] + "/QC_trimmed/final_multiqc_report.html"

rule trimmomatic:
    input: 
        R1 = config["data_dir"] + "/rawdata/P1658_{sample_full}_L002_R1_001_1.fq.gz",
        R2 = config["data_dir"] + "/rawdata/P1658_{sample_full}_L002_R2_001_2.fq.gz"
    output:
        trim_1P = config["proj_dir"] + "/trimmed/{sample_full}_R1_tr_1P.fq.gz",
        trim_1U = config["proj_dir"] + "/trimmed/{sample_full}_R1_tr_1U.fq.gz",
        trim_2P = config["proj_dir"] + "/trimmed/{sample_full}_R2_tr_2P.fq.gz",
        trim_2U = config["proj_dir"] + "/trimmed/{sample_full}_R2_tr_2U.fq.gz"
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
        config["proj_dir"] + "/trimmed/{sample_full}_{read}.fq.gz"
    output: 
        html = config["proj_dir"] + "/QC_trimmed/{sample_full}_{read}_fastqc.html",
        zip = config["proj_dir"] + "/QC_trimmed/{sample_full}_{read}_fastqc.zip"
    shell:
        """
        mkdir -p {config[proj_dir]}/QC_trimmed/
        fastqc {input} --outdir {config[proj_dir]}/QC_trimmed/
        """

rule multiqc_trimmed:
    input:
        expand(config["proj_dir"] + "/QC_trimmed/{sample_full}_{read}_fastqc.html", 
               sample_full=SAMPLES_FULL,
               read=["R1_tr_1P", "R2_tr_2P"])
    output:
        config["proj_dir"] + "/QC_trimmed/final_multiqc_report.html"
    params:
        report_name = "final_multiqc_report.html"
    shell:
        "multiqc {config[proj_dir]}/QC_trimmed/ -o {config[proj_dir]}/QC_trimmed/ -n {params.report_name}"
