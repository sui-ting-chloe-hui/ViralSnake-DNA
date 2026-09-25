# Define the reference path (use the decompressed version)
ref_path = config["ref_dir"] + config["ref_fasta"] 

rule all002:
    input:
        ref_path + ".fai",
        ref_path.replace(".fasta.gz", ".dict"),
        expand(ref_path + ".{ext}", ext=["bwt", "pac", "ann", "amb", "sa"])

rule samtools_faidx:
    input: ref_path
    output: ref_path + ".fai"
    shell: "gunzip -c {input} | samtools faidx - > {output}"

rule bwa_index:
    input: ref_path
    output:
        bwt = ref_path + ".bwt",
        pac = ref_path + ".pac",
        ann = ref_path + ".ann",
        amb = ref_path + ".amb", 
        sa = ref_path + ".sa"
    shell: "gunzip -c {input} | bwa index -p {input} -"

rule picard_dict:
    input: ref_path
    output: ref_path.replace(".fasta.gz", ".dict")
    shell:
        """
        gunzip -c {input} > temp_reference.fasta
        PicardCommandLine CreateSequenceDictionary \
          REFERENCE=temp_reference.fasta \
          OUTPUT={output}
        rm temp_reference.fasta
        """