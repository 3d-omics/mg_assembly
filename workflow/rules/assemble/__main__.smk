include: "__functions__.smk"
include: "megahit.smk"
include: "index.smk"
include: "bowtie2.smk"


rule assemble__all:
    """Run everything in the assemble module"""
    input:
        rules.assemble__megahit__all.input,
        rules.assemble__index__all.input,
        rules.assemble__bowtie2__all.input
