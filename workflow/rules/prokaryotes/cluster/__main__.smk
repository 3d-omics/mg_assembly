include: "__functions__.smk"
include: "concoct.smk"
include: "magscot.smk"
include: "maxbin2.smk"
include: "metabat2.smk"


rule prokaryotes__cluster:
    """Run the assemble module"""
    input:
        rules.prokaryotes__cluster__magscot__all.input,
