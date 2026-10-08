
## 002 - Reference indexing

# Define reference path once
ref_path = config["ref_dir"] + config["ref_fasta"]

rule all002:
    input:
        ref_path + ".fai",
        ref_path.rsplit(".", 1)[0] + ".dict",
        expand(ref_path + ".{ext}", ext=["bwt", "pac", "ann", "amb", "sa"])

rule samtools_faidx:
    input: ref_path
    output: ref_path + ".fai"
    shell: "samtools faidx {input}"

rule bwa_index:
    input: ref_path
    output:
        bwt = ref_path + ".bwt",
        pac = ref_path + ".pac", 
        ann = ref_path + ".ann",
        amb = ref_path + ".amb",
        sa = ref_path + ".sa"
    shell: "bwa index {input}"

rule picard_dict:
    input: ref_path
    output: ref_path.rsplit(".", 1)[0] + ".dict"
    shell:
        """
        PicardCommandLine CreateSequenceDictionary \
          REFERENCE={input} \
          OUTPUT={output}
        """