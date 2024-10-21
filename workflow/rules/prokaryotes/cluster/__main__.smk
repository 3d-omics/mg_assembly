include: "__functions__.smk"
include: "concoct.smk"
include: "magscot.smk"
include: "maxbin2.smk"
include: "metabat2.smk"


rule prokaryotes__cluster__all:
    """Run the assemble module"""
    input:
        rules.prokaryotes__cluster__concoct__all.input,
        rules.prokaryotes__cluster__maxbin2__all.input,
        rules.prokaryotes__cluster__metabat2__all.input,
        rules.prokaryotes__cluster__magscot__all.input,
