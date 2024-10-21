include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "coverm.smk"


rule prokaryotes__quantify__all:
    input:
        rules.prokaryotes__quantify__coverm__all.input,
