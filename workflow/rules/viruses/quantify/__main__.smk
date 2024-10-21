include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "coverm.smk"


rule viruses__quantify__all:
    input:
        rules.viruses__quantify__index__all.input,
        rules.viruses__quantify__bowtie2__all.input,
        rules.viruses__quantify__coverm__all.input,
