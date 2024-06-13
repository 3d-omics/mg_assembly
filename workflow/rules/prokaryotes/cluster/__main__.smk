include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "concoct.smk"
include: "drep.smk"
include: "magscot.smk"
include: "maxbin2.smk"
include: "metabat2.smk"


rule prokaryotes__cluster:
    """Run the assemble module"""
    input:
        rules.prokaryotes__cluster__drep.input,
