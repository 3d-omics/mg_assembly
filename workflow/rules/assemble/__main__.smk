include: "__functions__.smk"
include: "megahit.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "quast.smk"
include: "multiqc.smk"


rule assemble__all:
    """Run everything in the assemble module"""
    input:
        rules.assemble__megahit__all.input,
        rules.assemble__bowtie2__build__all.input,
        rules.assemble__bowtie2__map__all.input,
        rules.assemble__multiqc__all.input,
