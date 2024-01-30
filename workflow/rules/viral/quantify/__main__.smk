include: "__functions__.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "samtools.smk"


rule viral__quantify:
    input:
        rules.viral__quantify__coverm.input,
        rules.viral__quantify__samtools.input,
