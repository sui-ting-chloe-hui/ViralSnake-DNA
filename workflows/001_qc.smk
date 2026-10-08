## 001 - Initial quality control of raw fastq files
## Runs FastQC on every R1 and R2 file, then summarises with MultiQC.
## Raw file locations and naming come from raw_dir / r1_pattern / r2_pattern in the config.

import os
import re

# All raw read files for this project (R1 and R2 for every sample)
RAW_FASTQS = expand(RAW_R1, sample=SAMPLES) + expand(RAW_R2, sample=SAMPLES)

# FastQC names its output after the input file, minus .fastq.gz / .fq.gz
# e.g. CHIKV-1_R1_001.fastq.gz -> CHIKV-1_R1_001_fastqc.html
def fastqc_name(path):
    return re.sub(r"\.(fastq|fq)\.gz$", "", os.path.basename(path))

FASTQ_BY_NAME = {fastqc_name(fq): fq for fq in RAW_FASTQS}


rule all001:
    input:
        config["proj_dir"] + "/QC/raw_multiqc_report.html"


rule raw_fastqc:
    input:
        lambda wc: FASTQ_BY_NAME[wc.fq]
    output:
        html = config["proj_dir"] + "/QC/{fq}_fastqc.html"
    shell:
        """
        mkdir -p {config[proj_dir]}/QC/
        fastqc {input} --outdir {config[proj_dir]}/QC/
        """


rule multiqc:
    input:
        [config["proj_dir"] + "/QC/" + name + "_fastqc.html" for name in FASTQ_BY_NAME]
    output:
        config["proj_dir"] + "/QC/raw_multiqc_report.html"
    shell:
        """
        multiqc -f {config[proj_dir]}/QC/ -o $(dirname {output}) -n $(basename {output})
        """