include: "cluster/__main__.smk"
include: "annotate/__main__.smk"
include: "quantify/__main__.smk"
include: "multiqc.smk"

rule prokaryotes__all:
    input:
        rules.prokaryotes__cluster__all.input,
        rules.prokaryotes__annotate__all.input,
        rules.prokaryotes__quantify__all.input,
        rules.prokaryotes__multiqc__all.input,
