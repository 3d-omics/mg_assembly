include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"


rule quantify:
    input:
        rules.prokaryotes__quantify__coverm.input,
        rules.prokaryotes__quantify__samtools.input,
