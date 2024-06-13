include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"


rule viruses__quantify:
    input:
        rules.viruses__quantify__coverm.input,
        rules.viruses__quantify__samtools.input,
