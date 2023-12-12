include: "__functions__.smk"
include: "megahit.smk"
include: "bowtie2.smk"
include: "concoct.smk"
include: "metabat2.smk"
include: "maxbin2.smk"
include: "magscot.smk"
include: "drep.smk"


rule assemble__run:
    """Run all assembly rules (no evaluation)"""
    input:
        rules.assemble__drep.input,


rule assemble:
    """Run the assemble module"""
    input:
        rules.assemble__drep.input,
