include: "drep.smk"
include: "bowtie2.smk"


rule dereplicate_run:
    input:
        rules.dereplicate_bowtie2.input,
