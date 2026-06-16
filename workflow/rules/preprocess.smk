include: "preprocess/reads.smk"
include: "preprocess/hosts.smk"
include: "preprocess/fastp.smk"
include: "preprocess/bowtie2.smk"
include: "preprocess/clean.smk"
include: "preprocess/kraken2.smk"
include: "preprocess/bracken.smk"
include: "preprocess/singlem.smk"
include: "preprocess/nonpareil.smk"
include: "preprocess/multiqc.smk"


rule preprocess__all:
    input:
        rules.preprocess__reads__all.input,
        rules.preprocess__hosts__all.input,
        rules.preprocess__fastp__all.input,
        rules.preprocess__bowtie2__all.input,
        rules.preprocess__clean__all.input,
        rules.preprocess__kraken2__all.input,
        rules.preprocess__bracken__all.input,
        rules.preprocess__singlem__all.input,
        rules.preprocess__nonpareil__all.input,
        rules.preprocess__multiqc__all.input,
