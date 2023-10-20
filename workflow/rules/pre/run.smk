include: "fastp.smk"
include: "bowtie2.smk"


rule pre_run:
    input:
        rules.pre_bowtie2_extract_nonhost_all.input,
