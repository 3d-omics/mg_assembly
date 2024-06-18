include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "fastp.smk"
include: "kraken2.smk"
include: "nonpareil.smk"
include: "samtools.smk"
include: "singlem.smk"


rule preprocess:
    input:
        rules.preprocess__fastp.input,
        rules.preprocess__bowtie2.input,
        rules.preprocess__kraken2.input,
        rules.preprocess__nonpareil.output,
        rules.preprocess__singlem.input,
